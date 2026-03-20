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

  function _buildGenomeDropdown() {
    var current = _selectedGenome || _igvRef || '';
    var html = '<span style="font-size:12px;color:#888;font-weight:500;margin-right:4px">Reference:</span>';
    html += '<select id="__igv_genome_select__" style="font-size:12px;padding:4px 8px;max-width:220px;border:1px solid #ddd;border-radius:4px">';
    html += '<option value="' + (_igvRef || '') + '"' + (current === _igvRef ? ' selected' : '') + '>' + (_igvRef || 'none') + '</option>';
    KNOWN_GENOMES.forEach(function(g) {
      if (g.id !== _igvRef) {
        html += '<option value="' + g.id + '"' + (current === g.id ? ' selected' : '') + '>' + g.label + '</option>';
      }
    });
    html += '</select>';
    return html;
  }

  function _renderIgv(container, fileUrl, filename) {
    container.innerHTML = '<div id="__igv_div__" class="ap-loading">Loading...</div>';
    _loadIgvJs().then(function() {
      var div = document.getElementById('__igv_div__');
      if (!div) return;
      div.innerHTML = '';
      var activeRef = _selectedGenome || _igvRef;
      var opts = {};
      var knownIds = KNOWN_GENOMES.map(function(g) { return g.id; });
      if (knownIds.indexOf(activeRef) >= 0) {
        opts.genome = activeRef;
      } else {
        opts.reference = { fastaURL: '/file/' + encodeURIComponent(activeRef), indexed: false };
      }
      opts.tracks = [{ type: 'alignment', format: 'bam', url: fileUrl, name: filename }];
      igv.createBrowser(div, opts);
    }).catch(function(e) {
      container.innerHTML = '<div style="color:red;padding:16px;">IGV Error: ' + e.message + '</div>';
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
      _container = null;
      _metaCache = {};
      _currentFilename = '';
      delete window._bamPluginPaginate;
    }
  };
})();
