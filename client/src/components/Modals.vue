<template>
  <!-- TODO: should this file be incorporated to Home.vue? -->
  <div>
    <NewRequest/>
    <NewMeasurement/>
    <OpenProject/>
    <About/>
    <Login/>
    <GlobalSettings />
  </div>

</template>

<script>
  import GlobalSettings from './Modals/GlobalSettings.vue';
  import NewRequest from './Modals/NewRequest.vue';
  import NewMeasurement from './Modals/NewMeasurement.vue';
  import OpenProject from './Modals/OpenProject.vue';
  import About from './Modals/About.vue';
  import Login from './Modals/Login/Login.vue';

  export default {
    name: 'Modals',
    components: {
      GlobalSettings,
      About,
      OpenProject,
      NewMeasurement,
      NewRequest,
      Login,
    },
  };
</script>
