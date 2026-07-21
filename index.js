// AutoPipe Plugin: bam-viewer
// BAM alignment viewer with server-side pagination via samtools
// Uses /data/ endpoint for paginated data (Docker or local samtools)

(function() {
  var PAGE_SIZE = 100;
  var _container = null;
  var _metaCache = {};
  var _filterChrom = '';
  var _filterMinMapq = 0;
  var _showHeader = false;
  var _currentFilename = '';

  var FLAG_BITS = [
    { bit: 0x1, name: 'paired' },
    { bit: 0x2, name: 'proper' },
    { bit: 0x4, name: 'unmapped' },
    { bit: 0x8, name: 'mate_unmap' },
    { bit: 0x10, name: 'reverse' },
    { bit: 0x20, name: 'mate_rev' },
    { bit: 0x40, name: 'read1' },
    { bit: 0x80, name: 'read2' },
    { bit: 0x100, name: 'secondary' },
    { bit: 0x200, name: 'failQC' },
    { bit: 0x400, name: 'duplicate' },
    { bit: 0x800, name: 'supplementary' }
  ];

  var CIGAR_OPS = 'MIDNSHP=X';

  // ── IGV.js integration ──
  var KNOWN_GENOMES = [
    {id:'hg38', label:'Human (GRCh38/hg38)'},
    {id:'hg19', label:'Human (GRCh37/hg19)'},
    {id:'mm39', label:'Mouse (GRCm39/mm39)'},
    {id:'mm10', label:'Mouse (GRCm38/mm10)'},
    {id:'rn7',  label:'Rat (mRatBN7.2/rn7)'},
    {id:'rn6',  label:'Rat (Rnor_6.0/rn6)'},
    {id:'dm6',  label:'Fruit fly (BDGP6/dm6)'},
    {id:'ce11', label:'C. elegans (WBcel235/ce11)'},
    {id:'danRer11', label:'Zebrafish (GRCz11/danRer11)'},
    {id:'sacCer3',  label:'Yeast (sacCer3)'},
    {id:'tair10',   label:'Arabidopsis (TAIR10)'},
    {id:'galGal6',  label:'Chicken (GRCg6a/galGal6)'}
  ];
  var _igvRef = null;
  var _igvMode = 'data';
  var _selectedGenome = null;
  var _fileUrl = '';
  var _igvBrowser = null;

  // Width of the region IGV opens at. Alignment tracks stop rendering past
  // their 30kb visibility window; 2kb keeps individual reads legible as bands
  // rather than hairlines, which matters most on sparse WGS coverage.
  var IGV_WINDOW = 2000;

  function escapeHtml(str) {
    return str.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  function renderFlag(flag) {
    var f = parseInt(flag, 10);
    if (isNaN(f)) return flag;
    var html = '';
    for (var i = 0; i < FLAG_BITS.length; i++) {
      var on = (f & FLAG_BITS[i].bit) !== 0;
      if (on) {
        var cls = FLAG_BITS[i].name === 'unmapped' ? 'flag-unmapped' : 'flag-on';
        html += '<span class="flag-badge ' + cls + '">' + FLAG_BITS[i].name + '</span>';
      }
    }
    return html || '<span class="flag-badge flag-on">0</span>';
  }

  function renderCigar(cigar) {
    if (!cigar || cigar === '*') return '*';
    var html = '';
    var parts = cigar.match(/\d+[MIDNSHP=X]/g) || [];
    for (var i = 0; i < Math.min(parts.length, 20); i++) {
      var op = parts[i][parts[i].length - 1];
      html += '<span class="cigar-op cigar-' + op + '">' + parts[i] + '</span>';
    }
    if (parts.length > 20) html += '...';
    return html;
  }

  function renderMapq(q) {
    var v = parseInt(q, 10);
    if (isNaN(v)) return q;
    var cls = v >= 30 ? 'mapq-good' : v >= 10 ? 'mapq-ok' : 'mapq-bad';
    return '<span class="' + cls + '">' + v + '</span>';
  }

  async function fetchPage(name, page) {
    var resp = await fetch(
      '/data/' + encodeURIComponent(name) + '?page=' + page + '&page_size=' + PAGE_SIZE
    );
    return await resp.json();
  }

  function _fetchReference() {
    return fetch('/api/reference').then(function(r) { return r.json(); })
      .then(function(d) { _igvRef = d.reference || null; })
      .catch(function() { _igvRef = null; });
  }

  function _loadIgvJs() {
    return new Promise(function(resolve, reject) {
      if (window.igv) { resolve(); return; }
      var s = document.createElement('script');
      s.src = 'https://cdn.jsdelivr.net/npm/igv@3/dist/igv.min.js';
      s.onload = function() { resolve(); };
      s.onerror = function() { reject(new Error('Failed to load igv.js')); };
      document.head.appendChild(s);
    });
  }

  // The reference arrives either as a genome ID, a bare filename, or a full
  // path (show_results passes whatever the caller supplied). /file/ is keyed by
  // filename alone, so strip any directory part.
  function _refUrl(ref) {
    var base = String(ref).replace(/\\/g, '/').split('/').pop();
    return '/file/' + encodeURIComponent(base);
  }

  function _disposeIgvBrowser() {
    if (_igvBrowser) {
      // The host never calls destroy(), and igv.js keeps every browser it
      // creates in a module-level list. Without this, each tab switch leaks a
      // browser plus its listeners and caches.
      try { igv.removeBrowser(_igvBrowser); } catch (e) { /* already detached */ }
      _igvBrowser = null;
    }
  }

  // Probe for a sidecar index. igv.js otherwise guesses "<url>.bai" and fails
  // silently when the guess is wrong or the server never registered the index.
  // A 1-byte ranged GET is used rather than HEAD so this works on any server
  // that serves /file/.
  function _probeUrl(url) {
    return fetch(url, { headers: { Range: 'bytes=0-0' } })
      .then(function(r) { return r.ok ? url : null; })
      .catch(function() { return null; });
  }

  function _findIndex(fileUrl, exts) {
    var candidates = [];
    for (var i = 0; i < exts.length; i++) {
      candidates.push(fileUrl + '.' + exts[i]);                       // reads.bam.bai
      candidates.push(fileUrl.replace(/\.[^.\/]+$/, '.' + exts[i]));  // reads.bai
    }
    return candidates.reduce(function(chain, url) {
      return chain.then(function(found) { return found || _probeUrl(url); });
    }, Promise.resolve(null));
  }

  // ── Locus discovery ──
  //
  // Without an explicit locus igv.js opens at the whole first chromosome (or
  // the whole genome when the reference has several contigs), which is far
  // past the alignment track's visibility window — the track renders nothing.
  // So the opening view has to be anchored on the first actual alignment.
  //
  // Finding it must not go through /data/. That endpoint shells out to
  // samtools on the server, which the IGV tab otherwise never needs: igv.js
  // reads BGZF and the BAI itself. Reading the index here keeps the two tabs
  // independent, so IGV still works on a server without samtools or Docker.

  function _rangeFetch(url, start, end) {
    return fetch(url, { headers: { Range: 'bytes=' + start + '-' + end } })
      .then(function(r) {
        if (!r.ok) throw new Error('range request failed: ' + r.status);
        return r.arrayBuffer();
      });
  }

  function _inflateGzip(bytes) {
    var stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
    return new Response(stream).arrayBuffer();
  }

  // BGZF is a run of standalone gzip members. Each is decompressed on its own
  // so this does not rely on multi-member support in DecompressionStream.
  function _bgzfInflate(buffer) {
    var view = new DataView(buffer);
    var parts = [];
    var off = 0;

    function step() {
      if (off + 18 > buffer.byteLength) return Promise.resolve(parts);
      if (view.getUint8(off) !== 0x1f || view.getUint8(off + 1) !== 0x8b) {
        return Promise.resolve(parts);
      }
      // Locate the BC extra subfield, which carries the member length - 1.
      var xlen = view.getUint16(off + 10, true);
      var bsize = -1;
      var p = off + 12;
      var xend = p + xlen;
      while (p + 4 <= xend) {
        var slen = view.getUint16(p + 2, true);
        if (view.getUint8(p) === 66 && view.getUint8(p + 1) === 67) {
          bsize = view.getUint16(p + 4, true) + 1;
          break;
        }
        p += 4 + slen;
      }
      // A truncated trailing member just ends the walk — the caller only ever
      // needs the leading blocks of whatever range was fetched.
      if (bsize <= 0 || off + bsize > buffer.byteLength) return Promise.resolve(parts);

      var member = new Uint8Array(buffer, off, bsize);
      off += bsize;
      return _inflateGzip(member).then(function(out) {
        parts.push(new Uint8Array(out));
        return step();
      });
    }

    return step().then(function(list) {
      var total = 0;
      list.forEach(function(a) { total += a.length; });
      var out = new Uint8Array(total);
      var at = 0;
      list.forEach(function(a) { out.set(a, at); at += a.length; });
      return out;
    });
  }

  // @SQ names in BAM order — the BAI addresses references by this index.
  function _readBamRefNames(fileUrl) {
    return _rangeFetch(fileUrl, 0, 262143)
      .then(_bgzfInflate)
      .then(function(buf) {
        var dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
        if (dv.getUint32(0, true) !== 0x014D4142) throw new Error('not a BAM file');
        var off = 8 + dv.getInt32(4, true);          // magic + l_text + text
        var nRef = dv.getInt32(off, true);
        off += 4;
        var names = [];
        var decoder = new TextDecoder();
        for (var i = 0; i < nRef; i++) {
          var lName = dv.getInt32(off, true);
          off += 4;
          names.push(decoder.decode(new Uint8Array(buf.buffer, buf.byteOffset + off, lName - 1)));
          off += lName + 4;                          // name (NUL-terminated) + l_ref
        }
        return names;
      });
  }

  // First populated entry of the linear index → virtual offset of the first
  // alignment in the file.
  function _readBaiFirstOffset(indexUrl) {
    return fetch(indexUrl)
      .then(function(r) {
        if (!r.ok) throw new Error('index fetch failed: ' + r.status);
        return r.arrayBuffer();
      })
      .then(function(buf) {
        var dv = new DataView(buf);
        if (dv.getUint32(0, true) !== 0x01494142) throw new Error('not a BAI file');
        var off = 4;
        var nRef = dv.getInt32(off, true);
        off += 4;
        for (var r = 0; r < nRef; r++) {
          var nBin = dv.getInt32(off, true);
          off += 4;
          for (var b = 0; b < nBin; b++) {
            off += 4;                                // bin id
            var nChunk = dv.getInt32(off, true);
            off += 4 + nChunk * 16;                  // chunk begin/end pairs
          }
          var nIntv = dv.getInt32(off, true);
          off += 4;
          for (var i = 0; i < nIntv; i++) {
            var lo = dv.getUint32(off, true);
            var hi = dv.getUint32(off + 4, true);
            off += 8;
            if (lo !== 0 || hi !== 0) {
              // 64-bit virtual offset: high 48 bits address the BGZF block,
              // low 16 the position inside its uncompressed payload.
              return { refId: r, coffset: hi * 65536 + Math.floor(lo / 65536), uoffset: lo % 65536 };
            }
          }
        }
        throw new Error('index contains no aligned reads');
      });
  }

  function _readAlignmentAt(fileUrl, loc) {
    // A BGZF block holds at most 64KB; fetch several so the record cannot be
    // cut short by a block boundary.
    return _rangeFetch(fileUrl, loc.coffset, loc.coffset + 262143)
      .then(_bgzfInflate)
      .then(function(buf) {
        var dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
        return {
          refId: dv.getInt32(loc.uoffset + 4, true),
          pos: dv.getInt32(loc.uoffset + 8, true) + 1   // stored 0-based
        };
      });
  }

  function _locusAround(chrom, pos) {
    if (!chrom || chrom === '*' || !(pos >= 1)) return null;
    // Open just ahead of the first read so the window fills with data rather
    // than centring on the edge of the covered region.
    var start = Math.max(1, pos - Math.floor(IGV_WINDOW * 0.1));
    return chrom + ':' + start + '-' + (start + IGV_WINDOW);
  }

  function _resolveLocusFromIndex(fileUrl, indexUrl) {
    if (typeof DecompressionStream === 'undefined') return Promise.resolve(null);
    return Promise.all([_readBaiFirstOffset(indexUrl), _readBamRefNames(fileUrl)])
      .then(function(res) {
        var loc = res[0], names = res[1];
        return _readAlignmentAt(fileUrl, loc).then(function(rec) {
          var id = rec.refId >= 0 ? rec.refId : loc.refId;
          return _locusAround(names[id], rec.pos);
        });
      })
      .catch(function() { return null; });
  }

  // Fallback for BAMs with no usable index. Needs samtools on the server.
  function _resolveLocusFromData(filename) {
    return fetch('/data/' + encodeURIComponent(filename) + '?page=0&page_size=1')
      .then(function(r) { return r.json(); })
      .then(function(d) {
        var rec = d && d.rows && d.rows[0];
        if (!rec) return null;
        return _locusAround(rec[2], parseInt(rec[3], 10));   // RNAME, POS
      })
      .catch(function() { return null; });
  }

  function _resolveLocus(filename, fileUrl, indexUrl) {
    var viaIndex = indexUrl ? _resolveLocusFromIndex(fileUrl, indexUrl) : Promise.resolve(null);
    return viaIndex.then(function(locus) {
      return locus || _resolveLocusFromData(filename);
    });
  }

  function _buildGenomeDropdown() {
    var current = _selectedGenome || _igvRef || '';
    var refLabel = _igvRef ? String(_igvRef).replace(/\\/g, '/').split('/').pop() : '';
    var html = '<span style="font-size:12px;color:#888;font-weight:500;margin-right:4px">Reference:</span>';
    html += '<select id="__igv_genome_select__" style="font-size:12px;padding:4px 8px;max-width:220px;border:1px solid #ddd;border-radius:4px">';
    html += '<option value="' + escapeHtml(_igvRef || '') + '"' + (current === _igvRef ? ' selected' : '') + '>' + escapeHtml(refLabel || 'none') + '</option>';
    KNOWN_GENOMES.forEach(function(g) {
      if (g.id !== _igvRef) {
        html += '<option value="' + g.id + '"' + (current === g.id ? ' selected' : '') + '>' + g.label + '</option>';
      }
    });
    html += '</select>';
    return html;
  }

  function _renderIgv(container, fileUrl, filename) {
    _disposeIgvBrowser();
    container.innerHTML = '';
    var div = document.createElement('div');
    div.className = 'ap-loading';
    div.textContent = 'Loading...';
    container.appendChild(div);

    var activeRef = _selectedGenome || _igvRef;
    var knownIds = KNOWN_GENOMES.map(function(g) { return g.id; });
    var isKnownGenome = knownIds.indexOf(activeRef) >= 0;

    return Promise.all([
      _loadIgvJs(),
      _findIndex(fileUrl, ['bai', 'csi']),
      isKnownGenome ? Promise.resolve(null) : _findIndex(_refUrl(activeRef), ['fai'])
    ]).then(function(probes) {
      var trackIndex = probes[1], refIndex = probes[2];
      // Locus discovery needs the index URL, so it runs after the probes.
      return _resolveLocus(filename, fileUrl, trackIndex).then(function(locus) {
        return { locus: locus, trackIndex: trackIndex, refIndex: refIndex };
      });
    }).then(function(results) {
      var locus = results.locus, trackIndex = results.trackIndex, refIndex = results.refIndex;
      // The user may have switched tabs while the probes were in flight.
      if (!div.isConnected) return;
      div.textContent = '';
      div.className = '';

      var opts = {};
      if (isKnownGenome) {
        opts.genome = activeRef;
      } else {
        opts.reference = { fastaURL: _refUrl(activeRef) };
        if (refIndex) {
          // Indexed means igv.js range-reads the FASTA instead of pulling the
          // whole file into memory — the difference between a few KB and the
          // entire reference.
          opts.reference.indexURL = refIndex;
          opts.reference.indexed = true;
        } else {
          opts.reference.indexed = false;
        }
      }
      if (locus) opts.locus = locus;

      var track = { type: 'alignment', format: 'bam', url: fileUrl, name: filename };
      if (trackIndex) {
        track.indexURL = trackIndex;
      } else {
        // igv.js needs an index to random-access a BAM. Say so, rather than
        // leaving an empty track with no explanation.
        var note = document.createElement('div');
        note.style.cssText = 'padding:8px 12px;margin-bottom:8px;border-radius:4px;' +
          'background:#fff8e1;border:1px solid #ffe082;color:#795548;font-size:12px';
        note.textContent = 'No index (.bai) found next to this BAM. IGV cannot ' +
          'load alignments without one — run "samtools index" on the file.';
        container.insertBefore(note, div);
      }
      opts.tracks = [track];

      // Returned, not fire-and-forget: a rejected createBrowser used to become
      // an unhandled rejection and leave a blank pane with no explanation.
      return igv.createBrowser(div, opts).then(function(browser) {
        _igvBrowser = browser;
      });
    }).catch(function(e) {
      container.innerHTML = '<div style="color:red;padding:16px;">IGV Error: ' +
        escapeHtml(e && e.message ? e.message : String(e)) + '</div>';
    });
  }

  function renderTable(name, rows, total, page, refs) {
    var totalPages = Math.ceil(total / PAGE_SIZE) || 1;

    // Filter rows client-side by chrom and MAPQ
    var filtered = rows;
    if (_filterChrom || _filterMinMapq > 0) {
      filtered = rows.filter(function(rec) {
        // rec[2] = Chr, rec[4] = MAPQ
        if (_filterChrom && rec[2] !== _filterChrom) return false;
        if (_filterMinMapq > 0 && parseInt(rec[4], 10) < _filterMinMapq) return false;
        return true;
      });
    }

    // Extract unique chroms from refs or rows
    var chroms = [];
    if (refs && refs.length > 0) {
      chroms = refs.map(function(r) { return r.name; });
    } else {
      var seen = {};
      rows.forEach(function(rec) {
        if (rec[2] && !seen[rec[2]]) { seen[rec[2]] = true; chroms.push(rec[2]); }
      });
    }

    var html = '';

    // Summary
    html += '<div class="bam-summary">';
    html += '<span class="stat"><b>' + (total || 0).toLocaleString() + '</b> reads</span>';
    if (refs && refs.length > 0) {
      html += '<span class="stat"><b>' + refs.length + '</b> references</span>';
    }
    html += '</div>';

    // Header section
    var cached = _metaCache[name] || {};
    if (cached.header) {
      html += '<div class="bam-header-section">';
      html += '<div class="bam-header-toggle" id="bamHeaderToggle">' + (_showHeader ? '\u25BC' : '\u25B6') + ' BAM Header</div>';
      if (_showHeader) {
        html += '<div class="bam-header-content">' + escapeHtml(cached.header).replace(/\n/g, '<br>') + '</div>';
      }
      html += '</div>';
    }

    // Controls
    html += '<div class="bam-controls">';
    html += '<select id="bamChromFilter"><option value="">All chromosomes</option>';
    for (var ci = 0; ci < chroms.length; ci++) {
      html += '<option value="' + chroms[ci] + '"' + (chroms[ci] === _filterChrom ? ' selected' : '') + '>' + chroms[ci] + '</option>';
    }
    html += '</select>';
    html += '<input type="number" id="bamMapqFilter" placeholder="Min MAPQ" value="' + _filterMinMapq + '" min="0" max="255" style="width:90px;">';
    html += '</div>';

    // Table
    html += '<div class="bam-table-wrap" style="max-height:450px;overflow:auto;">';
    html += '<table class="bam-table"><thead><tr>';
    html += '<th>#</th><th>Read Name</th><th>Chr</th><th>Pos</th><th>MAPQ</th><th>CIGAR</th><th>Flags</th>';
    html += '</tr></thead><tbody>';

    var startIdx = page * PAGE_SIZE;
    for (var ri = 0; ri < filtered.length; ri++) {
      var rec = filtered[ri];
      // samtools view output: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
      html += '<tr>';
      html += '<td style="color:#aaa">' + (startIdx + ri + 1) + '</td>';
      html += '<td style="font-family:monospace;font-size:10px;">' + escapeHtml(rec[0] || '') + '</td>';
      html += '<td><span class="chr-badge">' + escapeHtml(rec[2] || '*') + '</span></td>';
      html += '<td>' + escapeHtml(rec[3] || '') + '</td>';
      html += '<td>' + renderMapq(rec[4] || '0') + '</td>';
      html += '<td>' + renderCigar(rec[5] || '*') + '</td>';
      html += '<td>' + renderFlag(rec[1] || '0') + '</td>';
      html += '</tr>';
    }
    html += '</tbody></table></div>';

    // Pagination
    if (totalPages > 1) {
      var safeName = name.replace(/'/g, "\\'");
      html += '<div class="bam-pagination">';
      html += '<button onclick="window._bamPluginPaginate(\'' + safeName + "'," + (page - 1) + ')"' +
        (page <= 0 ? ' disabled' : '') + '>&laquo; Prev</button>';
      html += '<span class="page-info">Page ' + (page + 1) + ' / ' + totalPages +
        ' (' + total.toLocaleString() + ' reads)</span>';
      html += '<button onclick="window._bamPluginPaginate(\'' + safeName + "'," + (page + 1) + ')"' +
        (page >= totalPages - 1 ? ' disabled' : '') + '>Next &raquo;</button>';
      html += '</div>';
    }

    return html;
  }

  function renderError(msg) {
    return '<div class="bam-plugin"><div class="bam-error">' +
      '<div style="font-size:48px;margin-bottom:16px;text-align:center">&#9888;</div>' +
      '<h3 style="text-align:center">Cannot Read BAM File</h3>' +
      '<p style="text-align:center">' + escapeHtml(msg) + '</p>' +
      '<p style="text-align:center;color:#888;font-size:12px">BAM is a binary format. To view data, the server needs <code>Docker</code> (auto-pulls samtools image) or <code>samtools</code> installed.</p>' +
      '</div></div>';
  }

  async function renderPage(name, page) {
    var target = _container;
    if (!target) return;

    // In tab mode, render into the content div
    var content = target.querySelector('#__plugin_content__');
    if (content) target = content;

    target.innerHTML = '<div class="ap-loading">Loading...</div>';

    var data = await fetchPage(name, page);
    if (data.error) {
      target.innerHTML = renderError(data.error);
      return;
    }

    // Cache metadata from first page
    if (page === 0) {
      _metaCache[name] = {
        header: data.header || '',
        refs: data.refs || [],
        col_headers: data.col_headers || []
      };
    }
    var cached = _metaCache[name] || {};

    var html = '<div class="bam-plugin">';
    html += renderTable(name, data.rows || [], data.total || 0, page, cached.refs);
    html += '</div>';
    target.innerHTML = html;

    // Bind events
    var ht = target.querySelector('#bamHeaderToggle');
    if (ht) ht.addEventListener('click', function() { _showHeader = !_showHeader; renderPage(name, page); });
    var cs = target.querySelector('#bamChromFilter');
    if (cs) cs.addEventListener('change', function() { _filterChrom = this.value; renderPage(name, page); });
    var mq = target.querySelector('#bamMapqFilter');
    if (mq) mq.addEventListener('change', function() { _filterMinMapq = parseInt(this.value, 10) || 0; renderPage(name, page); });
  }

  // Global pagination handler
  window._bamPluginPaginate = function(name, page) {
    if (page < 0) return;
    renderPage(name, page);
  };

  function _showView(container, fileUrl, filename) {
    // Every path through here replaces container.innerHTML, detaching any live
    // IGV browser — drop it before the DOM goes away.
    _disposeIgvBrowser();
    if (_igvRef) {
      var tabsHtml = '<div style="display:flex;gap:4px;margin-bottom:12px">';
      tabsHtml += '<button id="__tab_data__" style="padding:6px 16px;border:1px solid #ddd;border-radius:4px;cursor:pointer;font-size:13px;' + (_igvMode === 'data' ? 'background:#007bff;color:white;border-color:#007bff' : 'background:#f8f8f8') + '">Data</button>';
      tabsHtml += '<button id="__tab_igv__" style="padding:6px 16px;border:1px solid #ddd;border-radius:4px;cursor:pointer;font-size:13px;' + (_igvMode === 'igv' ? 'background:#007bff;color:white;border-color:#007bff' : 'background:#f8f8f8') + '">IGV</button>';
      tabsHtml += '</div>';
      if (_igvMode === 'igv') tabsHtml += _buildGenomeDropdown();
      container.innerHTML = tabsHtml + '<div id="__plugin_content__"></div>';

      container.querySelector('#__tab_data__').onclick = function() { _igvMode = 'data'; _showView(container, fileUrl, filename); };
      container.querySelector('#__tab_igv__').onclick = function() { _igvMode = 'igv'; _showView(container, fileUrl, filename); };
      var genomeSelect = container.querySelector('#__igv_genome_select__');
      if (genomeSelect) genomeSelect.onchange = function() { _selectedGenome = this.value; _showView(container, fileUrl, filename); };

      var content = container.querySelector('#__plugin_content__');
      if (_igvMode === 'igv') {
        _renderIgv(content, fileUrl, filename);
      } else {
        renderPage(filename, 0);
      }
    } else {
      renderPage(filename, 0);
    }
  }

  window.AutoPipePlugin = {
    render: function(container, fileUrl, filename) {
      // The host caches the plugin instance and only ever calls render(), so
      // this is the one reliable teardown point between files.
      _disposeIgvBrowser();
      _container = container;
      _currentFilename = filename;
      _fileUrl = fileUrl;
      _igvMode = 'data';
      _selectedGenome = null;
      _filterChrom = '';
      _filterMinMapq = 0;
      _showHeader = false;

      _container.innerHTML = '<div class="ap-loading">Loading...</div>';

      _fetchReference().then(function() {
        _showView(container, fileUrl, filename);
      });
    },
    destroy: function() {
      _disposeIgvBrowser();
      _container = null;
      _metaCache = {};
      _currentFilename = '';
      delete window._bamPluginPaginate;
    }
  };
})();
