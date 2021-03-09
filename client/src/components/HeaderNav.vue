<template>
  <div>
    <!-- Header-navbar -->
    <b-navbar type="dark" variant="dark">
      <b-navbar-nav>
        <b-navbar-brand>ePotku</b-navbar-brand>

        <!-- Site navigation links -->
        <b-nav-item-dropdown text="Navigation">
          <!-- Frontpage -->
          <b-dropdown-item to="/" :disabled="this.$route.path === '/'">Front page</b-dropdown-item>
          <!-- Main app -->
          <b-dropdown-item to="/app" :disabled="this.$route.path === '/app'">App</b-dropdown-item>
          <!-- TODO: Implement page and add :disabled="..." -->
          <b-dropdown-item>Open Data</b-dropdown-item>
          <!-- Admin page -->
          <b-dropdown-item
            v-if="currentUser.is_admin"
            to="/admin"
            :disabled="this.$route.path === '/admin'"
          >
            Admin
          </b-dropdown-item>
        </b-nav-item-dropdown>

        <b-nav-item-dropdown text="Help">
          <b-dropdown-item href="#">Manual</b-dropdown-item>
          <b-dropdown-item v-b-modal.about>About</b-dropdown-item>
        </b-nav-item-dropdown>

        <b-nav-item-dropdown text="File" v-if="this.$route.path === '/app'">
          <b-dropdown-item v-b-modal.newRequest>New Request</b-dropdown-item>
          <b-dropdown-item v-b-modal.newMeasurement>New Measurement</b-dropdown-item>
          <b-dropdown-divider></b-dropdown-divider>
          <b-dropdown-item v-b-modal.openProject>Open Request</b-dropdown-item>
        </b-nav-item-dropdown>

      </b-navbar-nav>

      <b-navbar-nav class="ml-auto" v-if="loading">
        <b-nav-item v-if="loading.tofe" @click="$emit('indexChange', 0)">
          <b-spinner small></b-spinner>  Loading...
        </b-nav-item>
        <b-nav-item v-else @click="$emit('indexChange', 0)">ToF-E histogram</b-nav-item>

        <b-nav-item v-if="loading.depth" @click="$emit('indexChange', 1)">
          <b-spinner small></b-spinner>  Loading...
        </b-nav-item>
        <b-nav-item v-else @click="$emit('indexChange', 1)">Depth profiles</b-nav-item>

        <b-nav-item v-if="loading.energy" @click="$emit('indexChange', 2)">
          <b-spinner small></b-spinner>  Loading...
        </b-nav-item>
        <b-nav-item v-else @click="$emit('indexChange', 2)">Energy spectra</b-nav-item>

        <b-nav-item v-if="loading.comp" @click="$emit('indexChange', 3)">
          <b-spinner small></b-spinner>  Loading...
        </b-nav-item>
        <b-nav-item v-else @click="$emit('indexChange', 3)">Element losses</b-nav-item>
      </b-navbar-nav>

      <!-- User info -->
      <b-navbar-nav class="ml-auto">
        <div v-if="isAuthenticated">
          <b-nav-item-dropdown :text="currentUser.email" right>
            <b-dropdown-item href="#">Profile (TODO)</b-dropdown-item>
            <b-dropdown-item @click="logout">Sign Out</b-dropdown-item>
          </b-nav-item-dropdown>
        </div>
        <div v-else>
          <b-nav-item v-b-modal.login>Log in</b-nav-item>
        </div>
      </b-navbar-nav>
    </b-navbar>
  </div>
</template>

<script>
  export default {
    name: 'HeaderNav',
    props: ['loading'],

    methods: {
      logout() {
        this.$store.dispatch('logout');
      },
    },

    computed: {
      currentUser() {
        return this.$store.state.user.attributes;
      },

      isAuthenticated() {
        return this.$store.getters.isAuthenticated;
      },
    },
  };
</script>
