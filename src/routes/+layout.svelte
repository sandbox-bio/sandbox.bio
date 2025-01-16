<script>
import { onMount } from "svelte";
import { page } from "$app/stores";
import {
	Styles,
	Navbar,
	Collapse,
	Nav,
	NavItem,
	NavLink,
	NavbarBrand,
	NavbarToggler,
	Container,
	Dropdown,
	DropdownToggle,
	DropdownMenu,
	DropdownItem,
	Modal,
	TabContent,
	TabPane,
	Tooltip,
	Toast,
	ToastHeader,
	ToastBody
} from "sveltestrap";
import Login from "$components/Login.svelte";
import LoginWithGoogle from "$components/LoginWithGoogle.svelte";
import { supabaseAnon } from "$src/utils";
import { user } from "$stores/user.js";
import { progress } from "$stores/progress.js";
import { env } from "$env/dynamic/public";

const playgrounds = [
	{ name: "Terminal", url: "/tutorials/playground", description: "Open ended" },
	{ name: "Jq", url: "/playgrounds/jq", description: "Wrangle JSON data" },
	{ name: "Awk", url: "/playgrounds/awk", description: "Wrangle tabular data" },
	{ name: "Grep", url: "/playgrounds/grep", description: "Search and filter" },
	{ name: "Sed", url: "/playgrounds/sed", description: "Search and replace" }
];

// -----------------------------------------------------------------------------
// State
// -----------------------------------------------------------------------------

let isNavbarOpen;
let loginModalOpen = false;
let toastOpen = false;
let toastToggle = () => (toastOpen = !toastOpen);

// Reactive state
export let data = {};
$: $progress = data.progress;
$: path = $page.url.pathname;

// -----------------------------------------------------------------------------
// Warn about losing progress if don't login
// -----------------------------------------------------------------------------

function remindLogin() {
	if ($user.email == null && path.startsWith("/tutorials")) toastToggle();
	setTimeout(remindLogin, 300000); // every 5 mins
}

// -----------------------------------------------------------------------------
// User auth
// -----------------------------------------------------------------------------

async function loginWithGoogle() {
	const redirectTo = $page.url.pathname;
	const result = await supabaseAnon.auth.signInWithOAuth({
		provider: "google",
		options: { redirectTo: `${$page.url.origin}/redirect?url=${redirectTo}` }
	});
	if (result.error) alert(result.error);
}

async function logout() {
	const data = await supabaseAnon.auth.signOut();
	if (data.error) console.error(data.error);
	else $user = {};
}

onMount(() => {
	setTimeout(remindLogin, 30000);
});
</script>

<svelte:head>
	<title>sandbox.bio</title>
	<script src="/v86/xterm.js"></script>

	<script async src="https://www.googletagmanager.com/gtag/js?id=G-DJ11EZ3RZ4"></script>
	<script>
	window.dataLayer = window.dataLayer || [];
	function gtag() {
		dataLayer.push(arguments);
	}
	gtag("js", new Date());
	gtag("config", "G-DJ11EZ3RZ4");
	</script>
</svelte:head>

<!-- Bootstrap CSS and icons -->
<Styles />

{#if env?.PUBLIC_USE_PRD_ASSETS}
	<div class="bg-primary-subtle text-primary text-center" style="font-size:10px">Development Mode</div>
{/if}

<!-- Navigation bar -->
<Navbar light container color="light" expand="md">
	<NavbarBrand href="/">&#129516; sandbox.bio</NavbarBrand>
	<NavbarToggler on:click={() => (isNavbarOpen = !isNavbarOpen)} />
	<Collapse isOpen={isNavbarOpen} navbar expand="md" on:update={(event) => (isNavbarOpen = event.detail.isOpen)}>
		<Nav class="ms-auto" navbar>
			<NavItem>
				<NavLink href="/tutorials" active={path.startsWith("/tutorials")}>Tutorials</NavLink>
			</NavItem>
			<Dropdown nav inNavbar>
				<DropdownToggle nav caret active={path.startsWith("/playgrounds")}>Playgrounds</DropdownToggle>
				<DropdownMenu end>
					{#each playgrounds as playground}
						<!--
							FIXME: Use "data-sveltekit-reload" so the browser reloads the page, not Svelte.
							This is to avoid cases where going to the Terminal playground, then to the other
							tools' playground makes the CodeMirror instances not responsive to typing, except
							for the backspace character.
						-->
						<DropdownItem href={playground.url} data-sveltekit-reload>
							<div class="d-flex">
								<span class="">{playground.name}</span>
								<span class="ps-2 text-muted opacity-75 small" style="padding-top:2px;">{playground.description}</span>
							</div>
						</DropdownItem>
					{/each}
				</DropdownMenu>
			</Dropdown>
			<NavItem>
				<NavLink href="/community" active={path.startsWith("/community")}>Community</NavLink>
			</NavItem>
			<NavItem>
				{#if $user?.email}
					<NavLink id="logout" on:click={logout}>Logout</NavLink>
					<Tooltip target="logout">
						<small>
							{$user.email}
						</small>
					</Tooltip>
				{:else}
					<NavLink on:click={() => (loginModalOpen = true)}>Log in</NavLink>
				{/if}
			</NavItem>
		</Nav>
	</Collapse>
</Navbar>

<!-- Toast Alert -->
{#if !env?.PUBLIC_USE_PRD_ASSETS}
	<div class="p-4 mb-4 me-3 position-fixed bottom-0 end-0" style="z-index: 15">
		<Toast autohide isOpen={toastOpen} header="" class="me-1">
			<ToastHeader toggle={toastToggle}>Note</ToastHeader>
			<ToastBody>Remember to log in to save your progress!</ToastBody>
		</Toast>
	</div>
{/if}

<!-- Login/Signup modal -->
<Modal body header="" toggle={() => (loginModalOpen = !loginModalOpen)} isOpen={loginModalOpen}>
	<TabContent>
		<!-- Login -->
		<TabPane tabId="login" active>
			<span class="h6" slot="tab">Log in</span>

			<!-- Login with Google -->
			<p class="mt-2 mb-2 small text-muted">Log in to save your progress:</p>
			<LoginWithGoogle direction="in" on:click={loginWithGoogle} />

			<!-- Or email/password -->
			<h6 class="mt-5">Or log in with your e-mail and password:</h6>
			<Login />
		</TabPane>

		<!-- Signup -->
		<TabPane tabId="signup">
			<span class="h6" slot="tab">Sign up</span>
			<p class="mt-2 mb-2 small text-muted">Create an account to save your progress:</p>
			<LoginWithGoogle direction="up" on:click={loginWithGoogle} />
		</TabPane>
	</TabContent>
</Modal>

<!-- Page Content -->
<Container class="mt-4">
	<slot />
</Container>

<!-- Footer (don't show on tutorials to avoid scrolling issues inside terminal) -->
{#if !$page.url.pathname.startsWith("/tutorials/")}
	<footer class="container pt-3 mt-5 mb-4 text-muted border-top">
		<div class="col-3">
			<h5>sandbox.bio</h5>
			<Nav vertical>
				<NavLink href="https://github.com/sandbox-bio/sandbox.bio/discussions" target="_blank" class="ps-0 py-1">Feedback</NavLink>
				<NavLink href="/about" class="ps-0 py-1">About</NavLink>
			</Nav>
		</div>
	</footer>

	<p class="mb-5" />
{/if}
