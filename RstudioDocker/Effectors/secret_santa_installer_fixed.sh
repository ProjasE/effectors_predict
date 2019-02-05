#! /bin/bash

# This file must be in the same folder as the tarballs.
# Path to the folder should be short (< 60 characters) #to do warings

### EDIT 2 LINES BELOW #########################################################
# SS_dep SS_out will be created in the directory with tarballs. Please replace 
# "SS_dep" with full path to the directory, e.g "/home/tools/SS_dep" and
# "SS_out" - accordingly.  
export SSDEP="/home/Effectors/SS_dep"
export SSOUT="/home/Effectors/SS_out"
################################################################################


export SSHDR="(SecretSanta Installer)"

function may_build {
  if [ ! -d "$1" ]; then
    if [ -e "$1" ]; then
      echo "$SSHDR Cannot create folder $(basename $1)"
      exit 2
    fi
    mkdir "$1"
  fi
}

may_build "$SSDEP"

# Archive files must be in the same folder as this script.
# I have some files with .tar.gz and some others with .tar.Z
find *.tar* -exec bash -c '\
  echo "$SSHDR Untar {}" && \
  tar -xzf "{}" -C "$SSDEP"' \
\;

# Unfortunately signalp versions 2 and 3 use sh while version 4 uses perl.
find "$SSDEP" -maxdepth 1 -type d -name "signalp*" -exec bash -c '\
  echo $SSHDR Updating {}/signalp && \
  sed -i \
    -e "s@SIGNALP=.*@SIGNALP='"'"'{}'"'"'@" \
    -e "s@ENV{SIGNALP} = .*@ENV{SIGNALP} = '"'"'{}'"'"'@" \
    -e "s@outputDir = .*@outputDir = '"'"'$SSOUT'"'"';@" "{}/signalp" && \
  mv "{}/signalp" "{}/signalp$(tr -cd 0-9 <<< "{}" | head -c1)"' \
\;

PERLDIR="$(which perl)"
GAWKDIR="$(which gawk)"

echo "$SSHDR Updating $SSDEP/targetp"
cd "$SSDEP"
TGPDIR="$(ls -d * | grep targetp)"
cd "$TGPDIR"
# caution, there is also TARGETPTMP later in targetp script!
sed -i \
  -e "s@TARGETP=.*@TARGETP='$SSDEP/$TGPDIR'@" \
  -e "s@^TMP=.*@TMP='$SSOUT'@" \
  -e "s@^PERL=.*@PERL='$PERLDIR'@" \
  -e "s@^AWK=.*@AWK='$GAWKDIR'@" targetp
cd ..

echo "$SSHDR Cloning WoLFPSort repository"
git clone https://github.com/fmaguire/WoLFPSort.git
cd WoLFPSort/bin
mv runWolfPsortSummary wolfpsort
cd ../..

echo "$SSHDR Updating tmhmm and tmhmmformat.pl"
cd tmhmm*/bin
sed -i "s@#!.*@#! $PERLDIR@" tmhmm
sed -i "s@#!.*@#! $PERLDIR@" tmhmmformat.pl
cd ../..

# -mindepth is not POSIX-compliant; better use find $PWD ! -path $PWD ...
# ~/.profile requires manual cleaning if this script is run several times.
echo "$SSHDR Updating ~/.profile"
EXTRADIR=$( \
  find $PWD -maxdepth 1 -mindepth 1 -type d \
  | sed -e 's@WoLFPSort@WoLFPSort/bin@' -e 's@\(tmhmm.*\)@\1/bin@' \
  | tr '\n' ':')
echo "export PATH="'"'"$EXTRADIR\$PATH"'"' >> ~/.profile
echo "IMPORTANT: now reload profile with . ~/.profile "
