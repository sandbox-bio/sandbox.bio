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
wrangler secret put SUPABASE_URL --env stg      # Supabase endpoint
wrangler secret put SUPABASE_API_KEY --env stg  # SECRET key available Supabase: Settings --> API
```
