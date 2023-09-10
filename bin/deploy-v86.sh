#!/bin/bash

set -u
aws --endpoint-url https://${CLOUDFLARE_ACCOUNT_ID}.r2.cloudflarestorage.com s3 sync --delete static/v86 s3://sandbox-bio/v86/
