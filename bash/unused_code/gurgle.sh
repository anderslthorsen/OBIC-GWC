#!/bin/bash

cd /data/OBIC/BIDS

# Find relative path for each file
for F in $(find -type f -name '*.nii.gz'); do
  # Parse out subject ID
  SUBJECT_ID=$(echo "$F" | cut -d"-" -f 2)

  # Move file, had sed replace -001 with "subject001" in path
  echo mv -v "$F" $(echo "$F" | sed "s/-$SUBJECT_ID/subject$SUBJECT_ID")

break
done
