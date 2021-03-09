<template>
  <!-- Sub-navbar for the main application -->
  <!-- TODO: Make this take less vertical spaces -->
  <b-navbar type="dark" variant="dark" style="">
    <b-navbar-nav>
      <b-nav-item-dropdown text="File">
        <b-dropdown-item v-b-modal.newProject>New Project</b-dropdown-item>
        <b-dropdown-item v-b-modal.newMeasurement>New Measurement</b-dropdown-item>
        <b-dropdown-divider />
        <b-dropdown-item v-b-modal.openProject>Open Project</b-dropdown-item>
      </b-nav-item-dropdown>
      <b-nav-item>View</b-nav-item>
      <b-nav-item>Tools</b-nav-item>
      <b-nav-item-dropdown text="Help">
        <b-dropdown-item href="#">Manual</b-dropdown-item>
        <b-dropdown-item v-b-modal.about>About</b-dropdown-item>
      </b-nav-item-dropdown>
    </b-navbar-nav>
  </b-navbar>
</template>
