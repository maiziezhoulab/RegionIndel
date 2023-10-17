# change permission
cd bin
chmod +x *.py
chmod -R 777 k8-0.2.4
tar -xzf SPAdes-3.13.0-Linux.tar.gz
rm SPAdes-3.13.0-Linux.tar.gz
cd ..

if ! [ -x "$(command -v samtools)" ];
then
    echo 'Error: samtools is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda samtools
else
    echo 'using existing samtools...'
fi

if ! [ -x "$(command -v minimap2)" ];
then
    echo 'Error: minimap2 is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda minimap2
else
    echo 'using existing minimap2...'
fi



echo 'You have installed AquilaSV dependencies and downloaded the source files successfully!'
 



