<!-- TODO: number formatting on focus change, e.g. 0 -> 0.000  -->
<!-- TODO: another modal for choosing the element from the periodic table? -->
<!-- TODO: prevent switching tabs if any field has an invalid value -->
<!-- TODO: red background on unfocused input fields with invalid values -->
<template>
  <b-modal :id="scope+'-settings'" :title="scope +' settings'" ok-title="Save" @ok="save"
           no-close-on-backdrop no-close-on-esc size="xl" scrollable>
    <b-tabs content-class="mt-3">
      <Measurement :scope="scope" :data="data.measurement"/>
      <Detector    :scope="scope" :data="data.detector"/>
      <Profile     :scope="scope" :data="data.profile"/>
      <Simulation  :scope="scope" :data="data.mcsimu" v-if="data.mcsimu"/>
    </b-tabs>
  </b-modal>
</template>

<script>
  /* eslint-disable no-param-reassign */
  import axios from 'axios';
  import Utils from '../../js/utils';
  import Measurement from './Tabs/Measurement.vue';
  import Detector from './Tabs/Detector.vue';
  import Simulation from './Tabs/Simulation.vue';
  import Profile from './Tabs/Profile.vue';

  export default {
    name: 'Settings',
    components: {
      Simulation,
      Detector,
      Measurement,
      Profile,
    },
    props: [
      'project',
      'scope',
    ],

    data() {
      return {
        data: {},
        okToClose: true,
      };
    },

    created() {
      if (this.scope === 'Request') {
        this.data = {
          measurement: {},
          detector: {},
          profile: {},
          mcsimu: {},
        };
      } else {
        this.data = {
          measurement: {},
          detector: {},
          profile: {},
        };
      }

      axios.get(`http://localhost:5000/settings/${this.scope}`)
        .then((res) => {
          Object.keys(res.data).forEach((setting) => {
            this.data[setting] = res.data[setting];
          });
        })
        .catch((error) => {
          // eslint-disable-next-line
          console.error(error);
        });
    },

    methods: {
      /**
       * Saves the user's changes to the settings.
       * TODO: prevent user from changing tabs while invalid inputs exist
       */
      save(modalEvent) {
        modalEvent.preventDefault();
        this.okToClose = true;

        Object.keys(this.data).forEach((key) => {
          if (!this.okToClose) return;
          this.iterate(this.data[key], key);
        });

        if (this.okToClose) {
          this.$nextTick(() => {
            // ensure backwards compatibility with desktop version by updating unused properties
            this.data.measurement.beam.charge = this.data.measurement.beam.charge_state;
            this.data.measurement.beam.ion = `${this.data.measurement.beam.isotope}${this.data.measurement.beam.element}`;
            this.$bvModal.hide(`${this.scope}-settings`);
          });
        }
      },

      /**
       * Updates each attribute in the object by calling iterate() until max depth
       * is reached (i.e. key no longer contains an object), then uses the HTML element's
       * ID to assign the user's input to the corresponding attribute.
       * TODO: v-model would probably be a lot simpler
       * TODO: cannot save added foils since iteration is through downloaded settings
       * TODO: saving foils wouldn't work anyway due to ID depth
       * (other settings are request-detector-setting, foils would be request-detector-foil-setting)
       *
       * @param obj Object to be iterated
       * @param tab Settings tab, used in getting the right <input> element
       */
      iterate(obj, tab) {
        Object.keys(obj).forEach((key) => {
          if (!this.okToClose) return; // don't bother finishing if even one setting is invalid
          if (typeof obj[key] === 'object') {
            this.iterate(obj[key], tab, this); // haven't reached max depth yet
          } else {
            const input = document.getElementById(`${this.scope}-${tab}-${key}`); // e.g. 'Request-measurement-name'
            if (input !== null) {
              if (key === 'modification_time') {
                obj[key] = (new Date()).toLocaleString('en-GB');
                obj[key.concat('_unix')] = Date.now(); // TODO: is unix time used for anything?
              } else {
                const inputIsInvalid = Utils.checkValidity(input, input.getAttribute('data-regex'));
                if (inputIsInvalid) {
                  this.okToClose = false;
                  return;
                }
                obj[key] = input.value; // assigns user's input to this.data
              }
            }
          }
        });
      },
    },
  };
</script>

<style>
  .column {
    width: 48%;
  }
</style>
