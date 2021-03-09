<template>
  <div id="composition-changes" class="maximize"></div>
</template>

<script>
  import Plotly from 'plotly.js-cartesian-dist';

  export default {
    name: 'ElementLosses',
    props: [
      'compositionChanges',
      'compLayout',
    ],
    methods: {
      draw() {
        if (this.compositionChanges.length) {
          Plotly.react(
            'composition-changes',
            this.compositionChanges,
            this.compLayout,
            { responsive: true },
          );
          document.getElementById('composition-changes').on('plotly_buttonclicked', (menu) => {
            if (menu.button.label === 'Save splits') { this.$emit('saveSplits'); }
            if (menu.button.name === 'reload') { this.$emit('reload'); }
          });
        }
      },
    },

    mounted() {
      this.draw();
    },
    watch: {
      compositionChanges() {
        this.draw();
      },
    },
  };
</script>
