<template>
  <b-modal id="openProject" hide-footer size="lg" title="Open an existing request ">
    <div class="d-block text-center">
      <Projects />
    </div>
  </b-modal>
</template>

<script>
  import Projects from '../Projects.vue';

  export default {
    name: 'OpenProject',
    components: {
      Projects,
    },
  };
</script>
