# Downloads the repository alongside the current repository.
# Creates softlinks to ensure responsiveness to edits to 'package'.

cd ..
if [ ! -d package ]
then
    git clone https://git.bioconductor.org/packages/simpleSingleCell package
fi
cd - 

ln -s ../package/vignettes/* .
