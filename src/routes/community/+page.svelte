<script>
import Alert from "$components/Alert.svelte";
import Listings from "$components/Listings.svelte";
import { config as IFB1 } from "$content/ifb-1/config.js";
import { config as IFB2 } from "$content/ifb-2/config.js";
</script>

<h3>Community Tutorials</h3>

<h5>French Institute of Bioinformatics</h5>
<img src="https://www.france-bioinformatique.fr/wp-content/uploads/logo-ifb-couleur.svg" alt="Logo of the IFB" width="300" />

<Listings title={null} items={[IFB1, IFB2]} showUnlisted={true} />

<!-- <div class="mt-5">
	<Alert>sdfsdf</Alert>
</div> -->
