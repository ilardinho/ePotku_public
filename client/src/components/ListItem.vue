<style scoped>
  .black_text a {
    padding: 8px 8px 8px 32px;
    text-decoration: none;
    color: #000000;
    display: block;
    transition: 0.3s;
  }

  .black_text a:hover {
    color: #515050;
  }

</style>

<template>
  <b-nav-item class="blac_text" @click="showModal(item.modal)">
    <b-icon v-bind:icon="item.icon" font-scale="1.5"></b-icon>{{item.text}}
  </b-nav-item>
</template>


<script>
  export default {
    name: 'ListItem',
    props: ['item'],
    methods: {
      showModal(id) {
        this.$bvModal.show(id);
      },
    },
  };
</script>
