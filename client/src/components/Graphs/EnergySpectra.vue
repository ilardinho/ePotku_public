<template>
  <div id="energy-spectra" class="maximize"></div>
</template>

<script>
  import Plotly from 'plotly.js-cartesian-dist';

  export default {
    name: 'EnergySpectra',
    props: [
      'energySpectra',
      'energyLayout',
    ],
    methods: {
      draw() {
        if (this.energySpectra.length) {
          Plotly.react(
            'energy-spectra',
            this.energySpectra,
            this.energyLayout,
            { responsive: true },
          );
        }
        document.getElementById('energy-spectra').on('plotly_buttonclicked', (menu) => {
          if (menu.button.name === 'reload') { this.$emit('reload'); }
        });
      },
    },

    mounted() {
      this.draw();
    },
    watch: {
      energySpectra() {
        this.draw();
      },
    },
  };
</script>
