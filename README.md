# sandbox.bio

## Infrastructure

|Environment|Domain|Access|Supabase DB|
|-|-|-|-|
|dev|[dev.sandbox.bio](https://dev.sandbox.bio)|[Only me](https://dash.teams.cloudflare.com/77294754f453e7c64b6100ddcde89b84/access/apps)|[dev](https://app.supabase.io/project/bqjvxpdzkembvixymfae/editor/table)|
|stg|[stg.sandbox.bio](https://stg.sandbox.bio)|[Testers](https://dash.teams.cloudflare.com/77294754f453e7c64b6100ddcde89b84/access/apps)|[stg](https://app.supabase.io/project/rrwfplicenewptmeeteq/editor/table)|
|prd|[[prd.]sandbox.bio](https://prd.sandbox.bio)|Public|[prd](https://app.supabase.io/project/vjmttfnyctkivaeljytg/editor/table)|


## Deploy

```bash
# Build Svelte App + Deploy Cloudflare Worker
npm run deploy-dev
```

## Setup

### Secrets

```bash
wrangler secret put SUPABASE_URL --env prd      # Supabase endpoint
wrangler secret put SUPABASE_API_KEY --env prd  # SECRET key available Supabase: Settings --> API

wrangler secret put SUPABASE_URL --env stg      # Supabase endpoint
wrangler secret put SUPABASE_API_KEY --env stg  # SECRET key available Supabase: Settings --> API

wrangler secret put SUPABASE_URL --env dev      # Supabase endpoint
wrangler secret put SUPABASE_API_KEY --env dev  # SECRET key available Supabase: Settings --> API
```
