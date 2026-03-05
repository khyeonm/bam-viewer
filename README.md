# bam-viewer

BAM alignment viewer with BGZF decompression, FLAG decoding, CIGAR visualization, and MAPQ filtering.

## Features

- BGZF block decompression via pako.js
- BAM binary header and alignment record parsing
- FLAG bitfield decoding with colored badges (paired, mapped, reverse, etc.)
- CIGAR string visualization with color-coded operations (M, I, D, S, H, N)
- MAPQ score display with quality indicator
- Reference sequence list from header (@SQ lines)
- Read group and program metadata display
- Chromosome and MAPQ filters
- Pagination for large files (100 reads per page, first 2000 reads max)

## Supported Extensions

- `.bam`

## Installation

Install from the **Plugins** tab in the AutoPipe desktop app.
