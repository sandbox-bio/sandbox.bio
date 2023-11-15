# sandbox.bio

[![Tests](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml/badge.svg)](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml)

Interactive bioinformatics command-line tutorials.

## Development

### Branches

- `main`: Development branch, pushing there runs tests
- `stg`: Staging branch, merge into this branch to deploy to stg.sandbox.bio
- `prd`: Production branch, merge into this branch to deploy to sandbox.bio

### Local dev

```bash
npm run dev
```

### Tests

- `npm run test` will launch all the tests in headless way
- `npx playwright test --ui` opens Playwright UI

### Deploy

```bash
# Generate static assets (see https://github.com/sandbox-bio/v86/blob/master/NOTES.md)
git clone https://github.com/sandbox-bio/v86.git && cd v86
make all
make build/xterm.js

# Generate .bin files
cd tools/docker/debian/
./generate.sh

# Upload .bin files
export CLOUDFLARE_ACCOUNT_ID=ID_GOES_HERE
./bin/deploy-v86.sh
```

---

## Infrastructure

### Subdomains

| Environment | Domain                                       | Access                                                                                    |
| ----------- | -------------------------------------------- | ----------------------------------------------------------------------------------------- |
| stg         | [stg.sandbox.bio](https://stg.sandbox.bio)   | [Testers](https://dash.teams.cloudflare.com/77294754f453e7c64b6100ddcde89b84/access/apps) |
| prd         | [sandbox.bio](https://sandbox.bio) | Public                                                                                    |

| Environment variable      | Description                  |
| ------------------------- | ---------------------------- |
| `PUBLIC_SUPABASE_URL`     | Supabase URL                 |
| `PUBLIC_SUPABASE_API_KEY` | Supabase database public key |
| `SUPABASE_API_KEY`        | Supabase database admin key  |

### Database

| Table | Description                                      | Access |
| ----- | ------------------------------------------------ | ------ |
| logs  | Log all calls to `sandbox.bio/*`                 | RLS    |
| pings | Analytics for tutorial progress                  | RLS    |
| state | Save tutorial progress                           | RLS    |

Append `_stg` to table names for dev/stg environments.
