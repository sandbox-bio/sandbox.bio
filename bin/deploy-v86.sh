#!/bin/bash
# FIXME: this script was only tested on MacOS (md5 and sed might have different flags)
set -u

ENV=${1?Usage: ./deploy.sh [dev|stg|prd]}
DIR_SOURCE=static/v86
DIR_DESTINATION=s3://sandbox-bio/v86/$ENV
FILE_DEBIAN_STATE=debian-state-base.bin.zst
URL_ENDPOINT=https://${CLOUDFLARE_ACCOUNT_ID}.r2.cloudflarestorage.com

[[ "$ENV" != "dev" ]] && [[ "$ENV" != "stg" ]] && [[ "$ENV" != "prd" ]] && echo "Unexpected environment '$ENV': expected dev, stg, or prd" && exit

# Add checksum to file name as a cache buster
CHECKSUM=$(md5 -q $DIR_SOURCE/$FILE_DEBIAN_STATE)
FILE_DEBIAN_STATE_CHECKSUM=debian-state-$CHECKSUM.bin.zst
cp $DIR_SOURCE/$FILE_DEBIAN_STATE $DIR_SOURCE/$FILE_DEBIAN_STATE_CHECKSUM
sed -i.bak 's/export const DEBIAN_STATE_ID =.*/export const DEBIAN_STATE_ID = "'$CHECKSUM'";/' src/config.js
rm src/config.js.bak

# Deploy if not dev
[[ "$ENV" == "dev" ]] && exit
AWS_MAX_ATTEMPTS=3 aws --endpoint-url "$URL_ENDPOINT" s3 sync --only-show-errors ${DIR_SOURCE} ${DIR_DESTINATION} --exclude "*" --include $FILE_DEBIAN_STATE_CHECKSUM
AWS_MAX_ATTEMPTS=3 aws --endpoint-url "$URL_ENDPOINT" s3 sync --only-show-errors --size-only --delete ${DIR_SOURCE}/debian-9p-rootfs-flat ${DIR_DESTINATION}/debian-9p-rootfs-flat/
