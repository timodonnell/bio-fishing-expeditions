# Interrupted first attempt

The initial scanner stopped on duplicate protein identifiers in shard 8 after 7,439,698 records in eight completed shards. These are a subset of the subsequent full run, not additional coverage. The exact executed source is `scripts/archive/scan_bgzf_initial.py`. Artifacts were moved here before a fresh run; each checkpoint records the relocation and unchanged content hashes. Outputs without a matching checkpoint JSON were never counted as completed. The final scanner gives every record a stable ordinal, preserving repeated names and sequences.
