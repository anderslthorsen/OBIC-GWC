#!/bin/bash -x
exit

cd /data/OBIC/BIDS

# Find relative path for each directory
#for F in sub-*; do
  # Parse out subject ID
  #SUBJECT_ID=$(echo "$F" | cut -d"-" -f 3)

  # Move file, had sed replace -001 with "subject001" in path
   #mv -vi "$F" $(echo "$F" | sed "s/-$SUBJECT_ID/subject$SUBJECT_ID/g")

#done

# Find relative path for each file
for F in $(find -type f -name '*.nii.gz'); do
  # Parse out subject ID
  SUBJECT_ID=$(basename "$F" | cut -d"-" -f 3)

  # Move file, had sed replace -001 with "subject001" in path
  mv -vi "$F" $(echo "$F" | sed "s/-$SUBJECT_ID/subject$SUBJECT_ID/g")

done
