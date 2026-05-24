v0.1:
  - Open .bam
  - Decode BGZF stream sequentially
  - Parse BAM header
  - Iterate records
  - Decode read name, flags, reference ID, position, CIGAR, sequence, quality, and aux tags
  - Output SAM-like text for validation against samtools view
