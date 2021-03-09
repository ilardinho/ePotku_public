<template>
  <div id="depth-profile" class="maximize"></div>
</template>

<script>
  import Plotly from 'plotly.js-cartesian-dist';

  export default {
    name: 'DepthProfiles',
    props: [
      'absoluteProfiles',
      'relativeProfiles',
      'depthLayout',
    ],
    methods: {
      draw() {
        Plotly.react(
          'depth-profile',
          this.absoluteProfiles,
          this.depthLayout,
          { responsive: true },
        );
        // swap between absolute and relative presentation mode on the graph
        document.getElementById('depth-profile').on('plotly_buttonclicked', (menu) => {
          if (menu.button.label === 'Relative') {
            this.depthLayout.legend.title.text = 'Units: Percentages';
            Plotly.react('depth-profile', this.relativeProfiles, this.depthLayout);
          }
          if (menu.button.label === 'Absolute') {
            this.depthLayout.legend.title.text = 'Units: Atoms/cm²';
            Plotly.react('depth-profile', this.absoluteProfiles, this.depthLayout);
          }
          if (menu.button.name === 'reload') { this.$emit('reload'); }
        });
      },
    },

    mounted() {
      this.draw();
    },

    watch: {
      absoluteProfiles() {
        this.draw();
      },
    },
  };
</script>
