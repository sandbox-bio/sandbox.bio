# sandbox.bio

## Deploy

```bash
# Build Svelte app
npm run build

# Build Svelte App + Deploy Cloudflare Worker
npm run deploy
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
