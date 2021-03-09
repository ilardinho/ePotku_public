<!-- TODO: disallow non-unique names for foils -->
<!-- TODO: checkbox for timing foil -->
<!-- TODO: compatibility with ProjectSettings save() function -->
<template>
  <div class="detector-row">
    <!-- foil name -->
    <input type="text" :value="foil.name">
    <!-- distance from target -->
    <span>{{foil.distance + ' mm'}}</span>
    <input v-on:keypress="checkInput($event, 'double')"
           :value="foil.distance" type="number" step="0.001">
    <!-- type of foil -->
    <select>
      <option>Circular</option>
      <option>Rectangular</option>
    </select>
    <!-- TODO: v-if rectangular -->
    <input type="number" step="0.001" v-on:keypress="checkInput($event, 'double')"
           :value="foil.diameter">
    <input type="number" step="0.001" v-on:keypress="checkInput($event, 'double')"
           :value="foil.transmission">
  </div>
</template>

<script>
  import Utils from '../../../js/utils';

  export default {
    name: 'Foil',
    props: ['foil', 'distance'],
    methods: {
      checkInput(keypress, regex) {
        Utils.checkInput(keypress, regex);
      },
    },
  };
</script>

<style>

</style>
