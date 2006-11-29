#! /bin/bash
#
# make the files to do appropriate regression testing.
#
PAIRCOMP=$1
VERSION=$2
for i in delta gcm otx
do
  $PAIRCOMP $i-sp.txt $i-lv.txt 20 .9 $i-$VERSION.cmp
done
