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
	DropdownItem
} from "sveltestrap";
import { URL_ASSETS } from "$src/config";

$: path = $page.url.pathname;

const playgrounds = [
	{ name: "Terminal", url: "/playgrounds/cli", description: "Open ended" },
	{ name: "Jq", url: "/playgrounds/jq", description: "Wrangle JSON data" },
	{ name: "Awk", url: "/playgrounds/awk", description: "Wrangle tabular data" },
	{ name: "Grep", url: "/playgrounds/grep", description: "Search and filter" },
	{ name: "Sed", url: "/playgrounds/sed", description: "Search and replace" }
];

// import Login from "$components/Login.svelte";

// // State
// let toastOpen = false;
// let toastToggle = () => toastOpen = !toastOpen;
// let loginModalOpen = false;
// let loginError = false;
// let loginSuccess = false;
// let loginBusy = false;
// let signupError = false;
// let signupSuccess = false;

// // -----------------------------------------------------------------------------
// // Warn about losing progress if don't login
// // -----------------------------------------------------------------------------

// function remindLogin() {
// 	if($user === null && path.startsWith("/tutorials"))
// 		toastToggle();
// 	setTimeout(remindLogin, 300_000);  // every 5 mins
// }
// setTimeout(remindLogin, 30_000);

// // -----------------------------------------------------------------------------
// // User auth
// // -----------------------------------------------------------------------------

// async function signup(credentials) {
// 	loginBusy = true;
// 	const data = await $supabase.auth.signUp(credentials);
// 	signupError = data.error?.message;
// 	if(data.error)
// 		loginBusy = false;
// 	else
// 		signupSuccess = "Account successfully created. Check your email for the verification link.";
// }

// async function login(credentials) {
// 	loginBusy = true;
// 	const data = await $supabase.auth.signIn(credentials);
// 	loginError = data.error?.message;
// 	if(data.error) {
// 		loginBusy = false;
// 	} else {
// 		loginError = false;
// 		loginModalOpen = false;
// 		$user = data.user;
// 		window.location.reload();
// 	}
// }

// async function logout() {
// 	const data = await $supabase.auth.signOut();
// 	console.error(data.error);
// 	if(data.error == null)
// 		$user = null;

// 	window.location.reload();
// }

onMount(async () => {
	// await envInit();
});

let isNavbarOpen;
</script>

<svelte:head>
	<title>sandbox.bio</title>
	<script src="{URL_ASSETS}/v86/xterm.js"></script>
</svelte:head>

<!-- Bootstrap CSS and icons -->
<Styles />

<Navbar light container color="light" expand="md">
	<NavbarBrand href="/">
		sandbox.bio
		<a href="/about">
			<span class="badge bg-danger rounded-pill beta ms-1 py-1">BETA</span>
		</a>
	</NavbarBrand>
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
						<DropdownItem href={playground.url}>
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
		</Nav>
	</Collapse>
</Navbar>

<!-- Toast Alert -->
<!-- <div class="p-4 mb-4 me-3 position-fixed bottom-0 end-0" style="z-index: 15">
	<Toast autohide isOpen={toastOpen} header="" class="me-1">
		<ToastHeader toggle={toastToggle}>Note</ToastHeader>
		<ToastBody>
			Remember to log in to save your progress!
		</ToastBody>
	</Toast>
</div> -->

<!-- Login/Signup modal -->
<!-- <Modal body header="" toggle={() => loginModalOpen = !loginModalOpen} isOpen={loginModalOpen}>
	<TabContent>
		<TabPane tabId="login" active>
			<span slot="tab"><h5>Log in</h5></span>
			<p class="mt-2 mb-2"><small>Log in to save your progress:</small></p>

			<Login type="login" error={loginError} success={loginSuccess} on:login={event => login(event.detail)} busy={loginBusy} />
		</TabPane>
		<TabPane tabId="signup" >
			<span slot="tab"><h5>Sign up</h5></span>
			<p class="mt-2 mb-2"><small>Create an account to save your progress:</small></p>

			<Login type="signup" error={signupError} success={signupSuccess} on:signup={event => signup(event.detail)} busy={loginBusy} />
		</TabPane>
	</TabContent>
</Modal> -->

<!-- Page Content -->
<Container class="mt-4">
	<slot />
</Container>

<!-- Footer -->
{#if !path.startsWith("/tutorials/") && path !== "/"}
	<footer class="container pt-3 mt-5 mb-4 text-muted border-top">
		<div class="col-3">
			<h5>sandbox.bio</h5>
			<Nav vertical>
				<NavLink href="https://github.com/sandbox-bio/sandbox.bio/discussions" target="_blank" class="ps-0 py-1">Feedback</NavLink>
				<NavLink href="/about" class="ps-0 py-1">About</NavLink>
			</Nav>
		</div>
	</footer>
{/if}

<style>
.beta {
	font-size: 10px;
	position: absolute;
}
</style>