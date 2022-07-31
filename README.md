# sandbox.bio

[![Deploy sandbox.bio](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy.yml/badge.svg)](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy.yml) [![Deploy stg.sandbox.bio](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-stg.yml/badge.svg)](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-stg.yml) [![Deploy dev.sandbox.bio](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-dev.yml/badge.svg)](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-dev.yml) [![Tests](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml/badge.svg)](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml)


## Infrastructure

### Domains

|Environment|Domain|Access|Supabase DB|
|-|-|-|-|
|dev|[dev.sandbox.bio](https://dev.sandbox.bio)|[Only me](https://dash.teams.cloudflare.com/77294754f453e7c64b6100ddcde89b84/access/apps)|[prd](https://app.supabase.io/project/vjmttfnyctkivaeljytg/editor/table)|
|stg|[stg.sandbox.bio](https://stg.sandbox.bio)|[Testers](https://dash.teams.cloudflare.com/77294754f453e7c64b6100ddcde89b84/access/apps)|[prd](https://app.supabase.io/project/vjmttfnyctkivaeljytg/editor/table)|
|prd|[[prd.]sandbox.bio](https://prd.sandbox.bio)|Public|[prd](https://app.supabase.io/project/vjmttfnyctkivaeljytg/editor/table)|


### Database

* https://app.supabase.io/project/vjmttfnyctkivaeljytg/editor/table

|Table|Description|Access|
|-|-|-|
|logs|Log all calls to sandbox.bio/*|RLS|
|pings|Analytics for tutorial progress|RLS|
|state|Save environment variables and tutorial progress|RLS|

Append `_stg` to table names for dev/stg environments.


## How to

### Create a new table

- [ ] Create table in dev
- [ ] In pgAdmin, right click on the table --> Scripts --> CREATE --> paste it in the UI for the new environment
- [ ] Enable Row-Level Security


### Deploy

```bash
# Build Svelte App + Deploy Cloudflare Worker
npm run deploy-dev
```

### Setup Secrets

```bash
wrangler secret put SUPABASE_URL --env prd      # Supabase endpoint
wrangler secret put SUPABASE_API_KEY --env prd  # SECRET key available Supabase: Settings --> API

wrangler secret put SUPABASE_URL --env stg      # Supabase endpoint
wrangler secret put SUPABASE_API_KEY --env stg  # SECRET key available Supabase: Settings --> API

wrangler secret put SUPABASE_URL --env dev      # Supabase endpoint
wrangler secret put SUPABASE_API_KEY --env dev  # SECRET key available Supabase: Settings --> API
```
