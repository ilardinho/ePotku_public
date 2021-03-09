<template>
  <b-modal ref="new-request-modal" id="newRequest" hide-footer title="Create a new Request">
    <!-- Here AddRequest vue is loaded into modal -->
    <AddRequest />
  </b-modal>
</template>

<script>
  import AddRequest from '../AddRequest.vue';

  export default {
    name: 'NewRequest',
    components: {
      AddRequest,
    },
  };
</script>
