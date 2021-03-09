<template>
  <!-- This is the center area of ePotku, this is where the graphs are supposed to be -->
  <b-tabs class="maximize" v-model="activeTab">

    <b-tab class="maximize" title="Tof-E histogram" lazy>
      <Heatmap v-if="tofeData.length"
               :id="'heatmap'"
               :ticks="ticks"
               :selections="elementSelections"
               :tofeData="tofeData"
               :compressionX="compressionX"
               :compressionY="compressionY"
               :xMin="xMin"
               :yMin="yMin"
               :transpose="transpose"
               v-on:getCuts="getCuts"
               v-on:reload="$refs['tofe-settings'].show()"/>
      <b-jumbotron v-else class="jumbotron">
        <template v-slot:header>Nothing here...</template>
        <template v-slot:lead>Click below to generate a ToF-E histogram.</template>
        <b-button v-b-modal.tofe-settings size="lg">ToF-E</b-button>
      </b-jumbotron>
    </b-tab>

    <b-tab title="Depth Profiles" class="maximize" lazy>
      <DepthProfiles v-if="absoluteProfiles.length"
                     :absoluteProfiles="absoluteProfiles"
                     :relativeProfiles="relativeProfiles"
                     :depthLayout="depthLayout"
                     v-on:reload="$refs['depth-settings'].show()"/>
      <b-jumbotron v-else class="jumbotron">
        <template v-slot:header>Nothing here...</template>
        <template v-slot:lead>Click below to generate a depth profiles graph.</template>
        <b-button v-b-modal.depth-settings size="lg">Depth profiles</b-button>
      </b-jumbotron>
    </b-tab>

    <b-tab title="Energy Spectra" class="maximize" lazy>
      <EnergySpectra v-if="energySpectra.length"
                     :energySpectra="energySpectra"
                     :energyLayout="energyLayout"
                     v-on:reload="$refs['energy-settings'].show()"/>
      <b-jumbotron v-else class="jumbotron">
        <template v-slot:header>Nothing here...</template>
        <template v-slot:lead>Click below to generate an energy spectra graph.</template>
        <b-button v-b-modal.energy-settings size="lg">Energy spectra</b-button>
      </b-jumbotron>
    </b-tab>

    <b-tab title="Element Losses" class="maximize" lazy>
      <ElementLosses v-if="compositionChanges.length"
                     :compositionChanges="compositionChanges"
                     :compLayout="compLayout"
                     v-on:reload="$refs['comp-settings'].show()"
                     v-on:saveSplits="saveSplits" />
      <b-jumbotron v-else class="jumbotron">
        <template v-slot:header>Nothing here...</template>
        <template v-slot:lead>Click below to generate an element losses graph.</template>
        <b-button v-b-modal.comp-settings size="lg">Element losses</b-button>
      </b-jumbotron>
    </b-tab>

    <!-- TODO: save values from modals to request/measurement settings -->
    <b-modal ref="depth-settings" id="depth-settings" @ok="getDepth" hide-header>
      <p v-if="!availableCuts.length">No available cuts</p>
      <div v-else class="depth-settings-content">
        <div class="column">
          <label class="settings-row"
                 v-for="cut in this.availableCuts" :key="cut">{{cut}}
            <input type="checkbox" :value="cut" v-model="selectedCuts">
          </label>
        </div>
        <div class="column">
          <!-- TODO: finish this bit -->
          <p>Settings info from request.PROFILE, potku2.ini, measurement.PROFILE,
            request.DETECTOR go here</p>
          <label class="settings-row">X axis units
            <select v-model="xUnits">
              <option>1e15 at./cm³</option>
              <option>nm</option>
              <!-- TODO: reference density option for nm scale, saved in measurement.profile -->
            </select>
          </label>
          <label class="settings-row">Systematic error:
            <input type="number" v-model="systematicError">
          </label>
          <label class="settings-row">Cross section
            <select v-model="crossSection">
              <option :value="1">Rutherford</option>
              <option :value="2">L'Ecuyer</option>
              <option :value="3">Andersen</option>
              <!-- TODO: reference density option for nm scale, saved in measurement.profile -->
            </select>
          </label>
          <label class="settings-row">Number of iterations:
            <input type="number" v-model="numOfIterations">
          </label>
        </div>
      </div>
    </b-modal>

    <b-modal ref="comp-settings" id="comp-settings"
             @ok="getCompositionChanges" hide-header size="sm">
      <p v-if="!availableCuts.length">No available cuts</p>
      <div v-else >
        <label class="settings-row"
               v-for="cut in this.availableCuts" :key="cut">{{cut}}
          <input type="checkbox" :value="cut" v-model="selectedCuts">
        </label>
        <label class="settings-row">Reference:
          <select v-model="compChangeSettings.referenceCut">
            <option disabled value="">Choose reference cut</option>
            <option v-for="cut in availableCuts" :key="cut">{{cut}}</option>
          </select>
        </label>
        <label class="settings-row">Splits:
          <input type="number" v-model="compChangeSettings.splits">
        </label>
        <label class="settings-row">Y-axis scale
          <select v-model="compChangeSettings.yAxisScale">
            <option>0 - max</option>
            <option>min - max</option>
          </select>
        </label>
      </div>
    </b-modal>

    <b-modal ref="energy-settings" id="energy-settings"
             @ok="getEnergySpectra" hide-header size="sm">
      <p v-if="!availableCuts.length">No available cuts</p>
      <div v-else >
        <label class="settings-row"
               v-for="cut in this.availableCuts" :key="cut">{{cut}}
          <input type="checkbox" :value="cut" v-model="selectedCuts">
        </label>
        <label class="settings-row">Histogram bin width:
          <input type="number" v-model="spectraBinWidth">
        </label>
      </div>
    </b-modal>

    <b-modal ref="tofe-settings" id="tofe-settings" @ok="getTofeData" size="sm" hide-header>
      <label class="settings-row">Compression x
        <input type="number" id="compression-x" v-model="compressionX">
      </label>
      <label class="settings-row">Compression y
        <input type="number" id="compression-y" v-model="compressionY">
      </label>
      <label class="settings-row">Transpose
        <input type="checkbox" v-model="transpose">
      </label>
    </b-modal>

  </b-tabs>
</template>

<script>
  import axios from 'axios';
  import Heatmap from './TofeHistogram.vue';
  import DepthProfiles from './Graphs/DepthProfiles.vue';
  import EnergySpectra from './Graphs/EnergySpectra.vue';
  import ElementLosses from './Graphs/ElementLosses.vue';

  export default {
    components: {
      ElementLosses,
      EnergySpectra,
      DepthProfiles,
      Heatmap,
    },

    props: [
      'loading',
      'activeTab',
    ],

    data() {
      return {
        // Tof-E data
        elementSelections: [],
        tofeData: [],
        ticks: {
          xVals: [],
          yVals: [],
          xText: [],
          yText: [],
        },
        xMin: 0,
        yMin: 0,
        transpose: true,

        // data for analysis graphs
        compositionChanges: [],
        energySpectra: [],
        absoluteProfiles: [],
        relativeProfiles: [],

        // ToF-E compression TODO: default values and changes should be saved in a settings file
        compressionX: 8,
        compressionY: 8,

        // composition changes settings
        compChangeSettings: {
          referenceCut: '',
          splits: 2,
          yAxisScale: '0 - max', // TODO: check how this is saved in the desktop version
        },

        // energy spectra settings
        spectraBinWidth: 0.025,

        // depth profile settings
        referenceDensity: 3.0, // measurement-specific .PROFILE file
        systematicError: 3.0, // widget_depth_profile.save
        xUnits: '1e15 at./cm³',
        crossSection: 1,
        crossSectionOptions: [
          { value: 1, text: 'Rutherford' },
          { value: 2, text: 'L´ecuyer' },
          { value: 3, text: 'Andersen' },
        ],
        numOfIterations: 4,

        // shared data
        selectedCuts: [],
        availableCuts: [],

        // energy spectra layout
        energyLayout: {
          title: 'Energy spectra',
          xaxis: { title: 'Energy (MeV)' },
          yaxis: { title: 'Yield (counts)' },
          legend: { font: { family: 'Lucida Console', size: 15 } },
          autosize: true,
          updatemenus: [{
            buttons: [{
              label: 'Reload',
              name: 'reload',
              method: 'skip',
            }],
            type: 'buttons',
            x: 1,
            xanchor: 'right',
            y: 1,
            yanchor: 'bottom',
            pad: {
              b: 10,
            },
          }],
        },

        // composition changes layout
        compLayout: {
          title: 'Composition changes', // TODO: display reference cut file
          xaxis: { title: 'Number of splits', autorange: true },
          yaxis: { title: 'Count', autorange: true },
          legend: { font: { family: 'Lucida Console', size: 15 } },
          autosize: true,
          updatemenus: [{
            buttons: [{
              label: 'Save splits',
              method: 'skip',
            }],
            type: 'buttons',
            pad: { b: 10 },
            x: 0,
            xanchor: 'left',
            y: 1,
            yanchor: 'bottom',
          }, {
            buttons: [{
              label: 'Reload',
              name: 'reload',
              method: 'skip',
            }],
            type: 'buttons',
            x: 1,
            xanchor: 'right',
            y: 1,
            yanchor: 'bottom',
            pad: {
              b: 10,
            },
          }],
        },

        // depth profile layout
        depthLayout: {
          autosize: true,
          title: 'Depth profile',
          legend: {
            title: {
              text: 'Units: Atoms/cm²',
              side: 'right',
              font: {
                family: 'Lucida Console',
                size: 15,
              },
            },
            font: {
              family: 'Lucida Console',
              size: 15,
            },
          },
          xaxis: {
            // TODO: set title based on depth units
            title: { text: 'Depth: 1e15 atoms/cm³' },
          },
          yaxis: {
            title: { text: 'Concentration' },
          },
          updatemenus: [{
            buttons: [{
              label: 'Absolute',
              method: 'skip',
            }, {
              label: 'Relative',
              method: 'skip',
            }],
            showactive: true,
            direction: 'right',
            type: 'buttons',
            pad: { b: 10 },
            x: 0,
            xanchor: 'left',
            y: 1,
            yanchor: 'bottom',
            font: { color: '#5072a8' },
          }, {
            buttons: [{
              label: 'Reload',
              name: 'reload',
              method: 'skip',
            }],
            type: 'buttons',
            x: 1,
            xanchor: 'right',
            y: 1,
            yanchor: 'bottom',
            pad: {
              b: 10,
            },
          }],
        },
      };
    },


    methods: {
      getCompositionChanges() {
        // TODO: notify user
        if (this.compChangeSettings.referenceCut === '') return;
        this.loading.comp = true;

        const compSettings = {
          referenceCut: this.compChangeSettings.referenceCut,
          yAxis: this.compChangeSettings.yAxisScale,
          splits: this.compChangeSettings.splits,
          cuts: this.selectedCuts,
        };

        axios.get('http://localhost:5000/element_losses', { params: compSettings })
          .then((res) => {
            if (Object.keys(res.data).length === 0) {
              this.loading.comp = false;
              return;
            }
            this.compGraph(res);
          })
          .catch((error) => {
            this.loading.comp = false;
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },


      compGraph(res) {
        // response data only has y values, x values is simply the number of splits
        const x = [...Array(res.data[0].data.length).keys()];

        // clear traces before adding new ones
        this.compositionChanges.length = 0;

        Object.keys(res.data).forEach((element) => {
          const change = res.data[element];
          this.compositionChanges.push({
            x,
            y: change.data,
            name: change.label,
            mode: 'lines',
            type: 'scatter',
            marker: { color: change.color },
            showlegend: true,
          });
        });
        this.loading.comp = false;
      },


      getEnergySpectra() {
        this.loading.energy = true;
        const energySettings = {
          binWidth: this.spectraBinWidth,
          cuts: this.selectedCuts,
        };
        axios.get('http://localhost:5000/load_energy_spectra', { params: energySettings })
          .then((res) => {
            this.energyGraph(res);
          })
          .catch((error) => {
            this.loading.energy = false;
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },


      energyGraph(res) {
        // clear traces before adding new ones
        this.energySpectra.length = 0;

        Object.keys(res.data).forEach((element) => {
          const spectra = res.data[element];
          this.energySpectra.push({
            x: spectra.x,
            y: spectra.y,
            // pad the start of each string to align the numbers on each row
            name: spectra.label,
            mode: 'lines',
            type: 'scatter',
            marker: { color: spectra.color },
            showlegend: true,
          });
        });
        this.loading.energy = false;
      },


      /**
       * Retrieves coordinates and draws depth profiles for each user-selected
       * element (this.selectedCuts)
       */
      getDepth() {
        this.loading.depth = true;
        // returns object, not array
        const depthSettings = {
          xUnits: this.xUnits,
          selectedCuts: this.selectedCuts,
          referenceDensity: this.referenceDensity,
          systematicError: this.systematicError,
        };
        axios.get('http://localhost:5000/load_depth', { params: depthSettings })
          .then((res) => {
            this.depthGraph(res);
          })
          .catch((error) => {
            this.loading.depth = false;
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },


      depthGraph(res) {
        // remove traces when recalculating concentrations
        this.relativeProfiles.length = 0;
        this.absoluteProfiles.length = 0;
        this.depthLayout.xaxis.title = `Depth: ${this.xUnits}`;

        // TODO: superscript for isotope (substring.sup())
        Object.keys(res.data.absolute).forEach((element) => {
          const profile = res.data.absolute[element];
          this.absoluteProfiles.push({
            x: profile.x,
            y: profile.y,
            // pad the start of each string to align the numbers on each row
            name: `${profile.label.padStart(4)} ${profile.absolute_concentration.toFixed(2).padStart(7)}`,
            // insert before cm2 to show depth in legend as well: atoms/1e153
            mode: 'lines',
            type: 'scatter',
            showlegend: true,
          });
        });

        Object.keys(res.data.relative).forEach((element) => {
          const profile = res.data.relative[element];
          this.relativeProfiles.push({
            x: profile.x,
            y: profile.y,
            // pad the start of each string to align the numbers on each row
            name: `${profile.label.padStart(4)} ${profile.percentage_concentration.toFixed(2).padStart(5)}% ± ${profile.error_margin.toFixed(2)}%`,
            mode: 'lines',
            type: 'scatter',
            showlegend: true,
          });
        });
        this.loading.depth = false;
      },


      /**
       * Retrieves the names of existing cut files
       */
      getCuts() {
        axios.get('http://localhost:5000/cuts')
          .then((res) => {
            this.availableCuts = res.data;
            this.selectedCuts = res.data;
          })
          .catch((error) => {
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },


      saveSplits() {
        if (this.compChangeSettings.referenceCut === '') return;

        const args = {
          referenceCut: this.compChangeSettings.referenceCut,
          splits: this.compChangeSettings.splits,
          cuts: this.selectedCuts,
        };

        axios.post('http://localhost:5000/element_losses', args)
          .then((res) => {
            // eslint-disable-next-line
            console.log(res);
          })
          .catch((error) => {
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },


      /**
       * Fetch data for ToF-E heatmap, including user selections
       */
      getTofeData() {
        // user input is saved as a string
        const compressionX = parseInt(this.compressionX, 10);
        const compressionY = parseInt(this.compressionY, 10);

        if (compressionX < 4 || compressionX > 20) return;
        if (compressionY < 4 || compressionY > 20) return;
        this.loading.tofe = true;
        // TODO: slider for compression values?
        axios.get('http://localhost:5000/selections', {})
          .then((res) => {
            /* ----------- BEGIN TOF-E ----------- */
            axios.get('http://localhost:5000/tofe', {
              params: { compressionX, compressionY },
              responseType: 'arraybuffer',
            })
              .then((tofe) => {
                const bytes = new Uint16Array(tofe.data);

                // first 8 bytes are metadata, set them to 0 after to prevent
                // them from showing up on the graph
                // eslint-disable-next-line prefer-destructuring
                this.xMin = bytes[0];
                bytes[0] = 0;
                // eslint-disable-next-line prefer-destructuring
                this.yMin = bytes[1];
                bytes[1] = 0;
                const xBins = bytes[2];
                bytes[2] = 0;
                const yBins = bytes[3];
                bytes[3] = 0;

                const tofeData = [];
                let byte = 0;

                // build ToF-E data array
                for (let i = 0; i < xBins; i++) {
                  const row = [];
                  for (let j = 0; j < yBins; j++) {
                    let num = bytes[byte++];
                    if (num > 0) {
                      if (num === 1) {
                        num = 1.5;
                      }
                      row.push(Math.log2(num));
                    } else {
                      row.push(0);
                    }
                  }

                  tofeData.push(row);
                }
                this.tofeData = tofeData;
                this.ticks.xVals.length = 0;
                this.ticks.yVals.length = 0;
                this.ticks.xText.length = 0;
                this.ticks.yText.length = 0;

                // scale range down according to histogram compression value
                // vals == X and Y values at which ticks are displayed on ToF-E graph
                // text == X and Y tick text at the corresponding values
                for (let i = 500; i <= 100000; i += 500) {
                  // ToF-E coordinate origin is (xMin, yMin), not (0,0) so we shift accordingly
                  this.ticks.xVals.push((i - this.xMin) / compressionX);
                  this.ticks.yVals.push((i - this.yMin) / compressionY);

                  this.ticks.xText.push(i);
                  this.ticks.yText.push(i);
                }
                /* ----------- END TOF-E ----------- */

                /* ----------- BEGIN SELECTIONS ----------- */
                this.elementSelections.length = 0;
                const sel = res.data;

                Object.keys(sel).forEach((key) => {
                  let name;
                  let element;

                  // TODO: need the isotope for scattered element (RBS) from server
                  if (sel[key].scatter) { // element type == RBS
                    element = sel[key].scatter;
                    // name = element + ('RBS').sup();
                    name = `${element}*`;
                  } else { // element type == ERD
                    element = sel[key].element;
                    name = ((sel[key].isotope).toString().sup() + element);
                  }

                  // also used in TofeHistogram.vue
                  const selection = {
                    x: [],
                    y: [],
                    type: 'scatter',
                    mode: 'markers+lines',
                    hoverinfo: 'none',
                    name,
                    element,
                    line: { color: sel[key].color },
                    marker: { size: 5, color: sel[key].color },
                    elementType: sel[key].element_type,
                    isotope: sel[key].isotope,
                    weightFactor: sel[key].weight_factor,
                  };

                  // scale coordinates down according to histogram compression value
                  // ToF-E coordinate origin is (xMin, yMin), not (0,0) so we shift accordingly
                  for (let i = 0; i < sel[key].points.length; i++) {
                    if (this.transpose) {
                      selection.x.push((sel[key].points[i][0] - this.xMin) / compressionX);
                      selection.y.push((sel[key].points[i][1] - this.yMin) / compressionY);
                    } else {
                      selection.y.push((sel[key].points[i][0] - this.xMin) / compressionX);
                      selection.x.push((sel[key].points[i][1] - this.yMin) / compressionY);
                    }
                  }

                  // add the first coords to the end so trace can be closed
                  selection.x.push(selection.x[0]);
                  selection.y.push(selection.y[0]);

                  // this.tofeHeatmap.push(selection);
                  this.elementSelections.push(selection);
                });
                /* ----------- END SELECTIONS ----------- */
              })
              .catch((error) => {
                this.loading.tofe = false;
                // eslint-disable-next-line no-console
                console.error(error);
              });

            this.loading.tofe = false;
          })
          .catch((error) => {
            this.loading.tofe = false;
            // eslint-disable-next-line
            console.log(error);
          });
      },
    },


    mounted() {
      // hides the b-tabs navigation bar
      const div = document.getElementsByClassName('tabs');
      div[0].children[0].style.display = 'none';

      // triggers caching on server
      this.getCuts();
    },
  };
</script>

<style>
  .depth-settings-content {
    display: flex;
    flex-direction: row;
    justify-content: space-between;
    /*flex-wrap: wrap;*/
  }


  .maximize {
    width: 100%;
    height: 100%
  }

  .tab-content {
    height: 100%;
  }


  .allow-overflow {
    overflow-x: hidden;
    height: 100%;
  }

  /* 42px is the tab selector height */
  .test {
    height: calc(100% - 42px);
  }

  .jumbotron {
    height: calc(100% - 42px);
    text-align: center;
  }
  .jumbotron template {
    vertical-align: middle;
  }
</style>
