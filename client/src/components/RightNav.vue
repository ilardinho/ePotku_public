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
  <b-nav class="bg-light sidebar">
    <!-- This is the part for collapsing the nav -->
    <b-navbar-toggle target="nav-collapse-right" :class="visible ? null : 'collapsed'"
      :aria-expanded="visible ? 'true' : 'false'"
      aria-controls="nav-collapse-right"
      @click="visible = !visible">
      <b-icon icon="justify" font-scale="1.5" class="sidebar" right></b-icon>

    </b-navbar-toggle>
    <b-collapse visible id="nav-collapse-right" is-nav v-model="visible">
    <!-- This is the part for collapsing the nav -->

      <ol class="nav flex-column">
        <b-nav class="black_text" vertical align="left">
          <ListItem v-for="item in data" v-bind:key="item.text" v-bind:item="item" />
        </b-nav>
      </ol>
    </b-collapse>
  </b-nav>
</template>

<script>
  // eslint-disable-next-line
  // import Axios from 'axios';
  import ListItem from './ListItem.vue';

  // const axios = require('axios').default;

  export default {
    name: 'RightNav',
    data() {
      return {
        visible: true,
        data: [{
                 icon: 'gear-wide',
                 text: 'Settings',
                 modal: 'projectSettings',
               },
               {
                 icon: 'screwdriver',
                 text: 'Measurement',
               },
               {
                 icon: 'gear-wide-connected',
                 text: 'Depth profile',
               },
               {
                 icon: 'hammer',
                 text: 'ToF calibration',
               },
               {
                 icon: 'brush',
                 text: 'Elemental selection',
               },
               {
                 icon: 'columns-gutters',
                 text: 'Save Cuts',
               },
               {
                 icon: 'graph-down',
                 text: 'Elemental Losses',
               },
               {
                 icon: 'graph-up',
                 text: 'Energy Spectrum',
               },
               {
                 icon: 'filter',
                 text: 'Depth Profile',
               },
        ],
      };
    },
    components: {
      ListItem,
    },
    /* axios tests */
    // created() {
    //
    //   const postTest = [{ icon: 'gear-wide', text: 'test' }];
    //   axios.post('http://papeniin1.student.it.jyu.fi/epotku_data', postTest)
    //     .then((response) => console.log(`post response ${response.data}`))
    //     .catch((error) => console.log(error));
    //
    //
    //   axios.get('http://papeniin1.student.it.jyu.fi/epotku_data.json')
    //     // axios.get('http://users.jyu.fi/~papeniin/ePotku/data.json') // CORS error
    //     .then((result) => { this.data = result.data; })
    //     // eslint-disable-next-line
    //     .catch((error) => console.log(error));
    // },
  };

</script>
