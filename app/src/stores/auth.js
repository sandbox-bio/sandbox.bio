import { readable, get } from "svelte/store";
import { createClient } from "@supabase/supabase-js";
import { config } from "./stores/config";
const $config = get(config);

const client = createClient($config.supabase.url, $config.supabase.publicKey);
export const supabase = readable(client);
