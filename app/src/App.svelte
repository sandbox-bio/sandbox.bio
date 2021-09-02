<script>
import { onMount } from "svelte";
import { Modal, TabContent, TabPane, Toast, ToastBody, ToastHeader, Icon } from "sveltestrap";
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";

import Home from "./routes/Home.svelte";
import Tutorial from "./tutorials/Tutorial.svelte";
import Terminal from "./terminal/Terminal.svelte";
import Login from "./components/Login.svelte";
import Listings from "./components/Listings.svelte";
import { config, supabase, user, progress, envInit } from "./stores/config";
import { tutorials } from "./stores/tutorials";

// Config
const intro = $config.playground;
const path = window.location.pathname;
const params = new URL(window.location).searchParams;

// State
let toastOpen = false;
let toastToggle = () => toastOpen = !toastOpen;
let loginModalOpen = false;
let loginError = false;
let loginSuccess = false;
let signupError = false;
let signupSuccess = false;

// -----------------------------------------------------------------------------
// Warn about losing progress if don't login
// -----------------------------------------------------------------------------

function remindLogin() {
	if($user === null && path.startsWith("/tutorials"))
		toastToggle();
	setTimeout(remindLogin, 60000);
}
setTimeout(remindLogin, 30000);


// -----------------------------------------------------------------------------
// User auth
// -----------------------------------------------------------------------------

async function signup(credentials) {
	const data = await $supabase.auth.signUp(credentials);
	signupError = data.error?.message;
	if(!data.error)
		signupSuccess = "Account successfully created. Check your email for the verification link.";
}

async function login(credentials) {
	const data = await $supabase.auth.signIn(credentials);
	loginError = data.error?.message;
	if(!data.error) {
		loginError = false;
		loginModalOpen = false;
		$user = data.user;
	}
}

async function logout() {
	const data = await $supabase.auth.signOut();
	console.error(data.error);
	if(data.error == null)
		$user = null;
}

onMount(async () => {
	await envInit();
})
</script>

<svelte:head>
	<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.4.1/font/bootstrap-icons.css">
</svelte:head>

<header class="d-flex flex-wrap justify-content-center py-3 border-bottom ms-4 me-4 mb-3">
	<a href="/" class="d-flex align-items-center mb-md-0 me-md-auto text-dark text-decoration-none">
		<span class="fs-4">🧬 sandbox.bio</span>&nbsp;
	</a>

	<ul class="nav nav-pills">
		<li class="nav-item dropdown">
			<!-- svelte-ignore a11y-invalid-attribute -->
			<a href="#" class="nav-link dropdown-toggle" class:active={path == "/tutorials"} id="navTutorials" role="button" data-bs-toggle="dropdown" aria-expanded="false">
				Tutorials
			</a>
			<ul class="dropdown-menu" aria-labelledby="navTutorials">
				<li><a class="dropdown-item" href="/tutorials">Browse all</a></li>
				<li><hr class="dropdown-divider"></li>
				{#each $tutorials as tutorial}
					<li>
						{#if $progress[tutorial.id]?.step == tutorial.steps.length - 1}
							<a class="dropdown-item text-success" href="/tutorials?id={tutorial.id}">
								<Icon name="check-circle-fill" />
								{tutorial.name}
							</a>
						{:else if $progress[tutorial.id]?.step > 0}
							<a class="dropdown-item text-primary" href="/tutorials?id={tutorial.id}&step={$progress[tutorial.id]?.step}">
								<Icon name="circle-half" />
								{tutorial.name}
							</a>
						{:else}
							<a class="dropdown-item" href="/tutorials?id={tutorial.id}">
								<Icon name="circle" />
								{tutorial.name}
							</a>
						{/if}
					</li>
				{/each}
			</ul>
		</li>
		<li class="nav-item">
			<a href="/playground" class="nav-link" class:active={path == "/playground"}>Playground</a>
		</li>
		<li class="nav-item">
			{#if $user == null}
				<button class="btn btn-link text-decoration-none" on:click={() => loginModalOpen = !loginModalOpen}>Log in</button>
			{:else}
				<div class="flex-shrink-0 dropdown ps-2">
					<!-- svelte-ignore a11y-invalid-attribute -->
					<a href="#" class="d-block link-dark text-decoration-none" id="dropdownUser" data-bs-toggle="dropdown" aria-expanded="false">
						<span class="text-primary" style="font-size: 1.7rem;">
							<Icon name="person-circle" />
						</span>
					</a>
					<ul class="dropdown-menu text-small shadow" aria-labelledby="dropdownUser" style="">
						<li><button class="dropdown-item disabled btn-sm" type="button">{$user.email}</button></li>
						<li><button class="dropdown-item" type="button" on:click={logout}>Log out</button></li>
					</ul>
				</div>
			{/if}
		</li>
	</ul>
</header>

<!-- Toast Alert -->
<div class="p-4 mb-4 me-3 position-fixed bottom-0 end-0" style="z-index: 15">
	<Toast autohide isOpen={toastOpen} header="" class="me-1">
		<ToastHeader toggle={toastToggle}>Note</ToastHeader>
		<ToastBody>
			Remember to log in to save your progress!
		</ToastBody>
	</Toast>
</div>

<!-- Main layout -->
<main role="main" class="ps-4 pe-4">
	{#if path == "/"}
		<Home />
	{:else if path.startsWith("/tutorials")}
		{#if params.get("id")}
			<Tutorial id={params.get("id")} step={+params.get("step") || 0} />
		{:else}
			<Listings items={$tutorials} />
		{/if}
	{:else if path.startsWith("/playground")}
		<div class="p-2" style="background-color:#000">
			<Terminal {intro} files={$tutorials[1].files} />
		</div>
	{/if}
</main>

<!-- Login/Signup modal -->
<Modal body header="" toggle={() => loginModalOpen = !loginModalOpen} isOpen={loginModalOpen}>
	<TabContent>
		<TabPane tabId="login" active>
			<span slot="tab"><h5>Log in</h5></span>
			<p class="mt-2 mb-2"><small>Log in to save your progress:</small></p>

			<Login type="login" error={loginError} success={loginSuccess} on:login={event => login(event.detail)} />
		</TabPane>
		<TabPane tabId="signup" >
			<span slot="tab"><h5>Sign up</h5></span>
			<p class="mt-2 mb-2"><small>Create an account to save your progress:</small></p>

			<Login type="signup" error={signupError} success={signupSuccess} on:signup={event => signup(event.detail)} />
		</TabPane>
	</TabContent>
</Modal>
