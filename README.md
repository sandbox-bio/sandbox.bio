# sandbox.bio

[![Deploy sandbox.bio](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy.yml/badge.svg)](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy.yml) [![Deploy stg.sandbox.bio](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-stg.yml/badge.svg)](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-stg.yml) [![Deploy dev.sandbox.bio](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-dev.yml/badge.svg)](https://github.com/sandbox-bio/sandbox.bio/actions/workflows/deploy-dev.yml) [![Tests](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml/badge.svg)](https://github.com/robertaboukhalil/sandbox.bio/actions/workflows/tests.yml)

Interactive bioinformatics command-line tutorials.

## Development

### Branches

- `dev`: Development branch, pushing there runs tests if changes to app/\* were made + auto deploys to dev.sandbox.bio
- `main`: Production branch; merge dev into main to deploy to stg.sandbox.bio

### Local dev

```bash
npm run dev
```

### Tests

- `npm run test` will launch all the tests in headless way
- `cypress open` opens Cypress so you can see the tests inside Chrome
- `npm run test -- --spec "tests/test_tutorials.js"` runs a single file

### Deploy

Deploys are done when committing to a branch, using Cloudflare Pages.

```bash
# Generate static assets (see https://github.com/sandbox-bio/v86/blob/master/NOTES.md)
git clone https://github.com/sandbox-bio/v86.git && cd v86
make all
make build/xterm.js

# Upload static assets
export CLOUDFLARE_ACCOUNT_ID=ID_GOES_HERE
./bin/deploy-v86.sh
```

### Example Quiz

```html
<script>
import Quiz from "$components/Quiz.svelte";
import Choice from "$components/QuizChoice.svelte";
</script>

<!-- The "id" is used to maintain state on page refresh. It just has to be unique within a step in a tutorial -->
<Quiz id="q1" choices={[
	{ valid: true, value: "Shell"},
	{ valid: false, value: "Terminal"},
]}>
	<span slot="prompt">
		Is Bash a shell or a terminal?
	</span>
</Quiz>

<Quiz id="q2" choices={[
	{ valid: true, value: "ls -s -h Data"},
	{ valid: true, value: "ls -sh Data"},
	{ valid: false, value: "ls -size -h Data"},
	{ valid: true, value: "ls --size -h Data"},
	{ valid: false, value: "ls --sizeh Data"},
	{ valid: false, value: "ls --size-h Data"},
	{ valid: true, value: "ls -h -s Data"},
	{ valid: true, value: "ls -hs Data"},
	{ valid: false, value: "ls -hsize Data"},
]}>
	<span slot="prompt">
		Among the following commands, which ones are correct?
	</span>
</Quiz>

<Quiz id="q3" choices={[
	{ valid: true, value: "yes"},
	{ valid: false, value: "no"},
]}>
	<span slot="prompt">
		Now, type the following command in your terminal and then press <kbd>Enter</kbd> key: `date`

		Does the terminal display the current date?
	</span>
</Quiz>
```

---

## Infrastructure

### Subdomains

| Environment | Domain                                       | Access                                                                                    |
| ----------- | -------------------------------------------- | ----------------------------------------------------------------------------------------- |
| dev         | [dev.sandbox.bio](https://dev.sandbox.bio)   | [Only me](https://dash.teams.cloudflare.com/77294754f453e7c64b6100ddcde89b84/access/apps) |
| stg         | [stg.sandbox.bio](https://stg.sandbox.bio)   | [Testers](https://dash.teams.cloudflare.com/77294754f453e7c64b6100ddcde89b84/access/apps) |
| prd         | [[prd.]sandbox.bio](https://prd.sandbox.bio) | Public                                                                                    |

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
| state | Save environment variables and tutorial progress | RLS    |

Append `_stg` to table names for dev/stg environments.
