#!/bin/bash

############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----                 GWAS dataset downloading                  -----    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

############    -----------------------------------------    ############
### ----------------------- About the script ------------------------ ###
############    -----------------------------------------    ############
# Author: Elizabeth Marquez-Gomez
# Date: 2022-September
# Version: 1
# Subversion: 0
# Using the file extension, this script will: determine download protocol to follow and renaming of the final file if possible.

#-----> Usage
# sh download.sh --url [path to url address] --location [local download path] --name [file name in source] --ext [file extension] --out [destination location for pre-proccessed GWAS]


############    -----------------------------------------    ############
### ------------------------------ Main ----------------------------- ###
############    -----------------------------------------    ############

if [ "$#" == "0" ]; then
echo "Error: No arguments" >&2; exit 1;
fi

while (( "$#" )); do
  case "$1" in
    -url)
      url=$2
      shift 2
      ;;
    -location)
      location=$2
      shift 2
      ;;
    -name)
      name=$2
      shift 2
      ;;
    -ext)
      ext=$2
      shift 2
      ;;
    -out)
      out=$2
      shift 2
      ;;
    -*|--*=|*) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

if [ $ext = 1 ]; then
  cd $location
  wget $url

  gunzip $name

  name=$(echo $name | sed 's/.gz//g')
  mv $name $out
  echo 'Download completed!'
fi
if [ $ext = 2 ]; then
  cd $location
  wget $url

  mv $name $out
  echo 'Download completed!'
fi
if [ $ext = 3 ]; then
  cd $location
  wget $url

  gzName=$(echo $name | sed 's/.bgz/.gz/g')
  mv $name $gzName
  gunzip $gzName

  name=$(echo $name | sed 's/.bgz//g')
  mv $name $out
  echo 'Download completed!'
fi
if [ $ext = 4 ]; then
  echo 'File needs special downloading'
  echo $url
  echo $out
fi
if [ $ext = 5 ]; then
  echo $name
  echo $out
  echo 'Download already done!'
fi
