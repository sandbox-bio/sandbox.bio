# sandbox.bio

[![Tests](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml/badge.svg)](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml)

Interactive bioinformatics command-line tutorials.

---

## Local development

### Branches

- `main`: Development branch, pushing there runs tests
- `stg`: Staging branch, merge into this branch to deploy to stg.sandbox.bio
- `prd`: Production branch, merge into this branch to deploy to sandbox.bio

### Environment setup

Define environment variables in `.env`:

```bash
# Supabase database URL
PUBLIC_SUPABASE_URL=...
# Supabase database public key
PUBLIC_SUPABASE_API_KEY=...
# Supabase database admin key
SUPABASE_API_KEY=...
```

Launch the web server:

```bash
npm install
npm run dev
```

### Tests

- `npm run test` will launch all the tests in headless way
- `npx playwright test --ui` opens Playwright UI

### Deploy

Generate Debian assets:

```bash
# Generate static assets (see https://github.com/sandbox-bio/v86/blob/master/NOTES.md)
git clone https://github.com/sandbox-bio/v86.git && cd v86
make all
make build/xterm.js

# Generate .bin files
cd tools/docker/debian/
./generate.sh
```

Deploy:

```bash
# Upload .bin files
export CLOUDFLARE_ACCOUNT_ID=ID_GOES_HERE

# Deploy main to stg
SOURCE=main
DEST=stg

# Deploy stg to prd
SOURCE=stg
DEST=prd

# Deploy assets and new debian image first (to avoid users caching a 404 for a non-existent image)
./bin/deploy-v86.sh $DEST
# Deploy code changes
git push origin --delete $DEST
git checkout -b $DEST --track origin/$SOURCE
git push origin $DEST
git checkout main
git branch -d $DEST
```

---

## Infrastructure

### Database

| Table | Description                      | Access |
| ----- | -------------------------------- | ------ |
| logs  | Log all calls to `sandbox.bio/*` | RLS    |
| pings | Analytics for tutorial progress  | RLS    |
| state | Save tutorial progress           | RLS    |

Append `_stg` to table names for dev/stg environments.

### Grafana

[Dashboard](https://sandboxbio.grafana.net)

- Import `grafana.json` into new Grafana instance: Depending on the setup, you might need to find/replace the postgres connection `uid` and replace it with the `uid` of the connection created.
- Export `grafana.json` to repo: Remove top-level `id`, `uid`, and `version` fields to avoid the import error `The dashboard has been changed by someone else`.
