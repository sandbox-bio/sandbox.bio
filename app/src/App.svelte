<script>
import { Modal, TabContent, TabPane } from "sveltestrap";
import { createClient } from "@supabase/supabase-js";
import "bootstrap/dist/js/bootstrap.bundle";
import "bootstrap/dist/css/bootstrap.min.css";

import Home from "./routes/Home.svelte";
import Tutorial from "./tutorials/Tutorial.svelte";
import Terminal from "./terminal/Terminal.svelte";
import Login from "./components/Login.svelte";
import Listings from "./components/Listings.svelte";
import { config } from "./stores/config";
import { tutorials } from "./stores/tutorials";

// Config
const intro = $config.playground;
const path = window.location.pathname;
const params = new URL(window.location).searchParams;
const supabase = createClient($config.supabase.url, $config.supabase.publicKey);

// State
let loginIsOpen = false;          // Whether "About" modal is showing or not
let userInfo = supabase.auth.user();  // Equals null if user isn't logged in
let loginError = false;
let loginSuccess = false;
let signupError = false;
let signupSuccess = false;

// -----------------------------------------------------------------------------
// User auth
// -----------------------------------------------------------------------------
async function signup(credentials) {
	const { user, session, error } = await supabase.auth.signUp(credentials);
	signupError = error?.message;
	if(!error)
		signupSuccess = "Account successfully created. Check your email for the verification link.";
}

async function login(credentials) {
	const { user, session, error } = await supabase.auth.signIn(credentials);
	loginError = error?.message;
	if(!error) {
		loginError = false;
		loginIsOpen = false;
		userInfo = user;
	}
}

async function logout() {
	const { error } = await supabase.auth.signOut();
	console.error(error);
	if(error == null)
		userInfo = null;
}
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
			<div class="input-group mb-3 nav-pills">
				<a role="button" class="btn btn-link text-decoration-none pe-1 nav-link" class:active={path == "/tutorials"} href="/tutorials">Tutorials</a>
				<button type="button" class="btn btn-link text-decoration-none dropdown-toggle dropdown-toggle-split ps-2 pe-3 nav-link" class:active={path == "/tutorials"} data-bs-toggle="dropdown" aria-expanded="false">
					<span class="visually-hidden">Toggle Dropdown</span>
				</button>
				<ul class="dropdown-menu">
					{#each $tutorials as tutorial}
						<li><a class="dropdown-item" href="/tutorials?id={tutorial.id}">{tutorial.name}</a></li>
					{/each}
				</ul>
			</div>
		</li>
		<li class="nav-item">
			<a href="/playground" class="nav-link" class:active={path == "/playground"}>Playground</a>
		</li>
		<li class="nav-item">
			{#if userInfo == null}
				<button class="btn btn-link text-decoration-none" on:click={() => loginIsOpen = !loginIsOpen}>Log in</button>
			{:else}
				<div class="flex-shrink-0 dropdown pt-1 ps-3">
					<a href="#" class="d-block link-dark text-decoration-none dropdown-toggle" id="dropdownUser2" data-bs-toggle="dropdown" aria-expanded="false">
						<img src="https://github.com/robertaboukhalil.png" alt="My profile" width="32" height="32" class="rounded-circle">
					</a>
					<ul class="dropdown-menu text-small shadow" aria-labelledby="dropdownUser2" style="">
						<li><button class="dropdown-item disabled btn-sm" type="button">{userInfo.email}</button></li>
						<li><button class="dropdown-item" type="button" on:click={logout}>Log out</button></li>
					</ul>
				</div>
			{/if}
		</li>
	</ul>
</header>

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
<Modal body header="" toggle={() => loginIsOpen = !loginIsOpen} isOpen={loginIsOpen}>
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

<!-- <button on:click={signup}>sdfsdf</button> -->
