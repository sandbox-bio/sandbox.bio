<script>
import { onMount } from "svelte";
import { Modal, TabContent, TabPane, Toast, ToastBody, ToastHeader, Icon } from "sveltestrap";
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";

import Home from "./routes/Home.svelte";
import Tutorial from "./tutorials/Tutorial.svelte";
import Login from "./components/Login.svelte";
import Listings from "./components/Listings.svelte";
import { supabase, user, progress, envInit } from "./stores/config";
import { tutorials } from "./stores/tutorials";

// Config
const path = window.location.pathname;
const params = new URL(window.location).searchParams;

// State
let toastOpen = false;
let toastToggle = () => toastOpen = !toastOpen;
let loginModalOpen = false;
let loginError = false;
let loginSuccess = false;
let loginBusy = false;
let signupError = false;
let signupSuccess = false;

// -----------------------------------------------------------------------------
// Warn about losing progress if don't login
// -----------------------------------------------------------------------------

function remindLogin() {
	if($user === null && path.startsWith("/tutorials"))
		toastToggle();
	setTimeout(remindLogin, 300_000);  // every 5 mins
}
setTimeout(remindLogin, 30_000);


// -----------------------------------------------------------------------------
// User auth
// -----------------------------------------------------------------------------

async function signup(credentials) {
	loginBusy = true;
	const data = await $supabase.auth.signUp(credentials);
	signupError = data.error?.message;
	if(data.error)
		loginBusy = false;
	else
		signupSuccess = "Account successfully created. Check your email for the verification link.";
}

async function login(credentials) {
	loginBusy = true;
	const data = await $supabase.auth.signIn(credentials);
	loginError = data.error?.message;
	if(data.error) {
		loginBusy = false;
	} else {
		loginError = false;
		loginModalOpen = false;
		$user = data.user;
		window.location.reload();
	}
}

async function logout() {
	const data = await $supabase.auth.signOut();
	console.error(data.error);
	if(data.error == null)
		$user = null;

	window.location.reload();
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
			<a href="#" class="nav-link dropdown-toggle" class:active={path == "/tutorials" && params.get("id") != "rosalind"} id="navTutorials" role="button" data-bs-toggle="dropdown" aria-expanded="false">
				Tutorials
			</a>
			<ul class="dropdown-menu" aria-labelledby="navTutorials">
				<li><a class="dropdown-item" href="/tutorials">Browse all</a></li>
				<li><hr class="dropdown-divider"></li>
				{#each $tutorials.filter(t => t.listed !== false && t.steps.length > 0) as tutorial, i}
					{#if tutorial.divider}
						<li><h6 class="dropdown-header mb-0 pb-1 { i > 0 ? "mt-2" : "" }">{tutorial.divider}</h6></li>
					{/if}
					<li>
						{#if $user !== null && $progress[tutorial.id]?.step == tutorial.steps.length - 1}
							<a class="dropdown-item text-success" href="/tutorials?id={tutorial.id}">
								<Icon name="check-circle-fill" />
								{tutorial.name}
							</a>
						{:else if $user !== null && $progress[tutorial.id]?.step > 0}
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
<main class="px-4">
	{#if path == "/"}
		<Home />
	{:else if path.startsWith("/tutorials")}
		{#if params.get("id")}
			<Tutorial id={params.get("id")} step={+params.get("step") || 0} />
		{:else}
			<Listings items={$tutorials} />
		{/if}
	{:else if path.startsWith("/playground")}
		<Tutorial id="playground" />
	{:else if path.startsWith("/rosalind")}
		<Tutorial id="rosalind" step={+params.get("step") || 0} />
	{/if}
</main>

<!-- Login/Signup modal -->
<Modal body header="" toggle={() => loginModalOpen = !loginModalOpen} isOpen={loginModalOpen}>
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
</Modal>
