# TODO:
# x Test 42bp, http labshare, s3/gcp URLs from 1K genomes
# - Test CRAM file subset --> mount EFS for CRAM ref genomes?

# ------------------------------------------------------------------------------
# Compile samtools/htslib on EC2 instance
# ------------------------------------------------------------------------------

ssh -i ~/Desktop/2021-08-17-test-aws.cer ec2-user@54.177.12.96
sudo yum install -y git zlib-devel bzip2-devel lzma liblzma-devel xz-devel libcurl-devel openssl-devel autoconf gcc ncurses-devel
sudo yum groupinstall "Development Tools"
git clone https://github.com/samtools/samtools.git
git clone https://github.com/samtools/htslib.git
cd htslib
git checkout 1.13
git submodule update --init --recursive
cd ../samtools
git checkout 1.13

autoheader
autoconf -Wno-syntax
./configure
make



# ------------------------------------------------------------------------------
# Bundle with Exodus
# TODO: look into installing musl library (doesn't seem to be picked up, even if set `LD_LIBRARY_PATH=/usr/local/musl/`)
# ------------------------------------------------------------------------------

# # Install musl
# cd ~
# git clone git://git.musl-libc.org/musl
# git checkout v1.2.2	
# cd musl
# ./configure
# make
# sudo make install

# Install Exodus
cd ~
sudo yum install -y pip
pip install --user exodus-bundler
exodus --tarball samtools/samtools | tar -zx
cd exodus/

# ------------------------------------------------------------------------------
# Create entrypoint
# ------------------------------------------------------------------------------

# Original code from https://docs.aws.amazon.com/lambda/latest/dg/runtimes-walkthrough.html
cat <<'EOF' > bootstrap
#!/bin/sh
set -euo pipefail
source $LAMBDA_TASK_ROOT/"$(echo $_HANDLER | cut -d. -f1).sh"
while true
do
  HEADERS="$(mktemp)"
  EVENT_DATA=$(curl -sS -LD "$HEADERS" -X GET "http://${AWS_LAMBDA_RUNTIME_API}/2018-06-01/runtime/invocation/next")
  REQUEST_ID=$(grep -Fi Lambda-Runtime-Aws-Request-Id "$HEADERS" | tr -d '[:space:]' | cut -d: -f2)

  URL_GCP="https://storage.googleapis.com/genomics-public-data/ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
  URL_AWS="https://1000genomes.s3.amazonaws.com/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
  URL_AWS_CRAM="https://1000genomes.s3.amazonaws.com/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam.cram"

  cd /tmp
  # /var/task/bin/samtools idxstats http://labshare.cshl.edu/shares/schatzlab/www-data/ribbon/DRR01684_merged.chr1.bam > /tmp/samtools.out 2>&1
  # /var/task/bin/samtools view -c http://labshare.cshl.edu/shares/schatzlab/www-data/ribbon/DRR01684_merged.chr1.bam 1:10e6-10.05e6 > /tmp/samtools.out 2>&1
  # /var/task/bin/samtools view -c "$URL_AWS_CRAM" 1:10e6-10.1e6 > /tmp/samtools.out 2>&1
  /var/task/bin/samtools view -H http://labshare.cshl.edu/shares/schatzlab/www-data/ribbon/DRR01684_merged.chr1.bam > /tmp/samtools.out 2>&1

  /var/task/bin/jq --arg RESPONSE "$(cat /tmp/samtools.out)" '{"statusCode":200,"body": $RESPONSE}' <<< '{}' > /tmp/response.json
  curl -X POST "http://${AWS_LAMBDA_RUNTIME_API}/2018-06-01/runtime/invocation/$REQUEST_ID/response" --data "@/tmp/response.json"
done
EOF

echo '' > function.sh

chmod 755 bootstrap function.sh
zip --symlinks -r9 ../samtools.zip *
aws lambda update-function-code --function-name samtools --zip-file fileb://../samtools.zip
time aws lambda invoke --function-name samtools --payload '{"text":"Hello"}' response.txt; cat response.txt





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Run on local Docker container: didn't work, probably b/c not using same OS!
# ------------------------------------------------------------------------------
# cd /src/samtools-lambda

# # Setup
# apt-get -y install musl musl-dev musl-tools
# pip install --user exodus-bundler
# #
# curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
# unzip awscliv2.zip
# ./aws/install
# rm awscliv2.zip
# rm -rf ./aws/


# # 
# exodus --tarball ./samtools/samtools | tar -zx  # /usr/local/bin/samtools
# cd exodus/

# # Deploy
# zip --symlinks -r9 ../samtools.zip *
# aws lambda update-function-code --function-name samtools --zip-file fileb://../samtools.zip

# # Launch
# time aws lambda invoke --function-name samtools --payload '{"text":"Hello"}' response.txt; cat response.txt
