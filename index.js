// AutoPipe Plugin: bam-viewer
// BAM alignment viewer with BGZF decompression via pako

(function() {
  var PAKO_CDN = 'https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js';
  var PAGE_SIZE = 100;
  var MAX_READS = 2000;
  var MAX_DECOMPRESS = 4 * 1024 * 1024; // 4MB

  var rootEl = null;
  var headerText = '';
  var refSeqs = [];
  var allReads = [];
  var filteredReads = [];
  var currentPage = 0;
  var filterChrom = '';
  var filterMinMapq = 0;
  var showHeader = false;

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

  function loadScript(url, cb) {
    if (window.pako) { cb(); return; }
    var s = document.createElement('script');
    s.src = url;
    s.onload = function() { cb(); };
    s.onerror = function() { cb(new Error('Failed to load pako')); };
    document.head.appendChild(s);
  }

  function decompressBGZF(buf) {
    var data = new Uint8Array(buf);
    var blocks = [];
    var totalSize = 0;
    var pos = 0;

    while (pos < data.length && totalSize < MAX_DECOMPRESS) {
      // Check BGZF magic
      if (data[pos] !== 31 || data[pos + 1] !== 139) break;
      // Get block size from BSIZE field (at offset 16 in extra field)
      var xlen = data[pos + 10] | (data[pos + 11] << 8);
      var bsize = -1;
      var extraPos = pos + 12;
      var extraEnd = extraPos + xlen;
      while (extraPos < extraEnd) {
        if (data[extraPos] === 66 && data[extraPos + 1] === 67) {
          bsize = (data[extraPos + 4] | (data[extraPos + 5] << 8)) + 1;
          break;
        }
        var slen = data[extraPos + 2] | (data[extraPos + 3] << 8);
        extraPos += 4 + slen;
      }
      if (bsize < 0) break;

      try {
        var block = pako.inflateRaw(data.subarray(pos + 18, pos + bsize - 8));
        blocks.push(block);
        totalSize += block.length;
      } catch(e) { break; }
      pos += bsize;
    }

    var result = new Uint8Array(totalSize);
    var offset = 0;
    for (var i = 0; i < blocks.length; i++) {
      result.set(blocks[i], offset);
      offset += blocks[i].length;
    }
    return result;
  }

  function parseBAM(decompressed) {
    var view = new DataView(decompressed.buffer, decompressed.byteOffset, decompressed.byteLength);
    var pos = 0;

    // Magic
    var magic = String.fromCharCode(decompressed[0], decompressed[1], decompressed[2], decompressed[3]);
    if (magic !== 'BAM\1') throw new Error('Not a valid BAM file');
    pos = 4;

    // Header
    var headerLen = view.getInt32(pos, true); pos += 4;
    var headerBytes = decompressed.subarray(pos, pos + headerLen);
    headerText = '';
    for (var i = 0; i < headerBytes.length; i++) headerText += String.fromCharCode(headerBytes[i]);
    pos += headerLen;

    // Reference sequences
    var nRef = view.getInt32(pos, true); pos += 4;
    refSeqs = [];
    for (var r = 0; r < nRef; r++) {
      var nameLen = view.getInt32(pos, true); pos += 4;
      var name = '';
      for (var ni = 0; ni < nameLen - 1; ni++) name += String.fromCharCode(decompressed[pos + ni]);
      pos += nameLen;
      var seqLen = view.getInt32(pos, true); pos += 4;
      refSeqs.push({ name: name, length: seqLen });
    }

    // Alignments
    allReads = [];
    while (pos + 4 < decompressed.length && allReads.length < MAX_READS) {
      if (pos + 36 > decompressed.length) break;
      var blockSize = view.getInt32(pos, true); pos += 4;
      if (blockSize <= 0 || pos + blockSize > decompressed.length) break;
      var blockStart = pos;

      var refID = view.getInt32(pos, true); pos += 4;
      var posn = view.getInt32(pos, true); pos += 4;
      var binMqNl = view.getUint32(pos, true); pos += 4;
      var nameLen2 = binMqNl & 0xff;
      var mapq = (binMqNl >> 8) & 0xff;
      var flagNc = view.getUint32(pos, true); pos += 4;
      var nCigarOp = flagNc & 0xffff;
      var flag = (flagNc >> 16) & 0xffff;
      var seqLen2 = view.getInt32(pos, true); pos += 4;
      pos += 4; // next refID
      pos += 4; // next pos
      pos += 4; // tlen

      var readName = '';
      for (var qi = 0; qi < nameLen2 - 1; qi++) readName += String.fromCharCode(decompressed[pos + qi]);
      pos += nameLen2;

      // CIGAR
      var cigarStr = '';
      for (var ci = 0; ci < nCigarOp; ci++) {
        var cigarVal = view.getUint32(pos, true); pos += 4;
        var opLen = cigarVal >> 4;
        var opCode = cigarVal & 0xf;
        cigarStr += opLen + (opCode < CIGAR_OPS.length ? CIGAR_OPS[opCode] : '?');
      }

      pos = blockStart + blockSize;

      allReads.push({
        name: readName,
        flag: flag,
        chrom: refID >= 0 && refID < refSeqs.length ? refSeqs[refID].name : '*',
        pos: posn + 1,
        mapq: mapq,
        cigar: cigarStr || '*'
      });
    }
  }

  function renderFlag(flag) {
    var html = '';
    for (var i = 0; i < FLAG_BITS.length; i++) {
      var on = (flag & FLAG_BITS[i].bit) !== 0;
      var cls = on ? (FLAG_BITS[i].name === 'unmapped' ? 'flag-unmapped' : 'flag-on') : 'flag-off';
      if (on) html += '<span class="flag-badge ' + cls + '">' + FLAG_BITS[i].name + '</span>';
    }
    return html || '<span class="flag-badge flag-on">0</span>';
  }

  function renderCigar(cigar) {
    if (cigar === '*') return '*';
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
    var cls = q >= 30 ? 'mapq-good' : q >= 10 ? 'mapq-ok' : 'mapq-bad';
    return '<span class="' + cls + '">' + q + '</span>';
  }

  function applyFilter() {
    filteredReads = allReads.filter(function(r) {
      if (filterChrom && r.chrom !== filterChrom) return false;
      if (r.mapq < filterMinMapq) return false;
      return true;
    });
    currentPage = 0;
  }

  function render() {
    if (!rootEl) return;
    var chroms = [];
    var seen = {};
    for (var i = 0; i < allReads.length; i++) {
      if (!seen[allReads[i].chrom]) { seen[allReads[i].chrom] = true; chroms.push(allReads[i].chrom); }
    }

    var totalPages = Math.max(1, Math.ceil(filteredReads.length / PAGE_SIZE));
    if (currentPage >= totalPages) currentPage = totalPages - 1;
    var startIdx = currentPage * PAGE_SIZE;
    var pageReads = filteredReads.slice(startIdx, startIdx + PAGE_SIZE);

    var html = '<div class="bam-plugin">';

    // Summary
    html += '<div class="bam-summary">';
    html += '<span class="stat"><b>' + filteredReads.length.toLocaleString() + '</b> reads</span>';
    html += '<span class="stat"><b>' + refSeqs.length + '</b> references</span>';
    if (allReads.length >= MAX_READS) html += '<span class="stat" style="color:#c62828">(first ' + MAX_READS + ' reads)</span>';
    html += '</div>';

    // Header section
    if (headerText) {
      html += '<div class="bam-header-section">';
      html += '<div class="bam-header-toggle" id="bamHeaderToggle">' + (showHeader ? '\u25BC' : '\u25B6') + ' BAM Header</div>';
      if (showHeader) {
        html += '<div class="bam-header-content">' + headerText.replace(/</g, '&lt;').replace(/\n/g, '<br>') + '</div>';
      }
      html += '</div>';
    }

    // Controls
    html += '<div class="bam-controls">';
    html += '<select id="bamChromFilter"><option value="">All chromosomes</option>';
    for (var ci = 0; ci < chroms.length; ci++) {
      html += '<option value="' + chroms[ci] + '"' + (chroms[ci] === filterChrom ? ' selected' : '') + '>' + chroms[ci] + '</option>';
    }
    html += '</select>';
    html += '<input type="number" id="bamMapqFilter" placeholder="Min MAPQ" value="' + filterMinMapq + '" min="0" max="255" style="width:90px;">';
    html += '</div>';

    // Table
    html += '<div class="bam-table-wrap" style="max-height:450px;overflow:auto;">';
    html += '<table class="bam-table"><thead><tr>';
    html += '<th>#</th><th>Read Name</th><th>Chr</th><th>Pos</th><th>MAPQ</th><th>CIGAR</th><th>Flags</th>';
    html += '</tr></thead><tbody>';

    for (var ri = 0; ri < pageReads.length; ri++) {
      var r = pageReads[ri];
      html += '<tr>';
      html += '<td style="color:#aaa">' + (startIdx + ri + 1) + '</td>';
      html += '<td style="font-family:monospace;font-size:10px;">' + r.name + '</td>';
      html += '<td><span class="chr-badge">' + r.chrom + '</span></td>';
      html += '<td>' + r.pos.toLocaleString() + '</td>';
      html += '<td>' + renderMapq(r.mapq) + '</td>';
      html += '<td>' + renderCigar(r.cigar) + '</td>';
      html += '<td>' + renderFlag(r.flag) + '</td>';
      html += '</tr>';
    }
    html += '</tbody></table></div>';

    // Pagination
    if (totalPages > 1) {
      html += '<div class="bam-pagination">';
      html += '<button data-page="prev">&laquo; Prev</button>';
      var startP = Math.max(0, currentPage - 3);
      var endP = Math.min(totalPages, startP + 7);
      if (startP > 0) html += '<button data-page="0">1</button><span>...</span>';
      for (var p = startP; p < endP; p++) {
        html += '<button data-page="' + p + '"' + (p === currentPage ? ' class="current"' : '') + '>' + (p + 1) + '</button>';
      }
      if (endP < totalPages) html += '<span>...</span><button data-page="' + (totalPages - 1) + '">' + totalPages + '</button>';
      html += '<button data-page="next">Next &raquo;</button>';
      html += '<span class="page-info">Page ' + (currentPage + 1) + ' of ' + totalPages + '</span>';
      html += '</div>';
    }

    html += '</div>';
    rootEl.innerHTML = html;

    // Events
    var ht = rootEl.querySelector('#bamHeaderToggle');
    if (ht) ht.addEventListener('click', function() { showHeader = !showHeader; render(); });
    var cs = rootEl.querySelector('#bamChromFilter');
    if (cs) cs.addEventListener('change', function() { filterChrom = this.value; applyFilter(); render(); });
    var mq = rootEl.querySelector('#bamMapqFilter');
    if (mq) mq.addEventListener('change', function() { filterMinMapq = parseInt(this.value, 10) || 0; applyFilter(); render(); });
    var pbs = rootEl.querySelectorAll('.bam-pagination button');
    for (var bi = 0; bi < pbs.length; bi++) {
      pbs[bi].addEventListener('click', function() {
        var pg = this.getAttribute('data-page');
        if (pg === 'prev') { if (currentPage > 0) currentPage--; }
        else if (pg === 'next') { if (currentPage < totalPages - 1) currentPage++; }
        else { currentPage = parseInt(pg, 10); }
        render();
      });
    }
  }

  window.AutoPipePlugin = {
    render: function(container, fileUrl, filename) {
      rootEl = container;
      rootEl.innerHTML = '<div class="bam-loading">Loading ' + filename + '...</div>';
      headerText = ''; refSeqs = []; allReads = []; filteredReads = [];
      currentPage = 0; filterChrom = ''; filterMinMapq = 0; showHeader = false;

      loadScript(PAKO_CDN, function(err) {
        if (err) { rootEl.innerHTML = '<div class="bam-error">Failed to load pako library.</div>'; return; }

        fetch(fileUrl)
          .then(function(resp) { return resp.arrayBuffer(); })
          .then(function(buf) {
            try {
              var decompressed = decompressBGZF(buf);
              parseBAM(decompressed);
              filteredReads = allReads.slice();
              render();
            } catch(e) {
              rootEl.innerHTML = '<div class="bam-error">Error parsing BAM: ' + e.message + '</div>';
            }
          })
          .catch(function(err) {
            rootEl.innerHTML = '<div class="bam-error">Error loading file: ' + err.message + '</div>';
          });
      });
    },
    destroy: function() { allReads = []; filteredReads = []; refSeqs = []; rootEl = null; }
  };
})();
