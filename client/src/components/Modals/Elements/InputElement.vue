<!-- TODO: Split this into separate componenets (InputElementText, InputElementButton, etc.) or
           replace completely with basic HTML input elements. -->
<!-- TODO: fix alignment for two Inputelements on the same row -->
<!-- TODO: add checkbox support -->
<!-- TODO: it would be better to replace this way of adding input elements with something cleaner-->
<template>
  <div class="settings-row">
    <!-- Clickable label -->
    <label v-if="input" :for="id" class="settings-label">{{title}}:</label>
    <!-- Not clickable -->
    <label v-else class="settings-label">{{title}}:</label>

    <!-- Button used for opening an element selection modal -->
    <button v-if="button" :id="id" class="settings-input" :value="button">{{button}}</button>
    <!-- Modification time, should not be editable -->
    <span v-else-if="time" :id="id" class="settings-input">{{time}}</span>
    <!-- Text area for custom descriptions -->
    <textarea v-else-if="textarea" :id="id" class="settings-input" :value="textarea"></textarea>
    <!-- Dropdown selections -->
    <select v-else-if="options" :id="id" :name="id" class="settings-input">
      <option selected="selected">{{selected}}</option>
      <option v-for="option in unselected" :key="option">{{option}}</option>
    </select>
    <!-- Disabled option -->
    <input v-else-if="disabled" class="settings-input" :value="input" disabled>
    <!-- The most common element, content is always a number (integer, float or e notation) -->
    <input v-else :id="id" class="settings-input" :value="input" :data-regex="regex"
           v-on:keypress="checkInput" v-on:change="checkValidity" :type="type" :step="step">
  </div>
</template>

<script>
  /* eslint-disable vue/script-indent */
  import Utils from '../../../js/utils';

  export default {
    name: 'InputElement',
    props: [
      'disabled',
      'title',
      'id',
      'input',
      'time',
      'textarea',
      'selected',
      'button',
      'options',
      'regex',
    ],
    methods: {
      checkInput(keypress) {
        Utils.checkInput(keypress, this.regex);
      },

      checkValidity(e) {
        Utils.checkValidity(e.target, this.regex);
      },

    },
    computed: {
      unselected() { // remove selected option from the list of options
        return this.options.filter((option) => option !== this.selected);
      },

      step() {
        if (this.regex === 'double') { return 0.001; }
        return 1;
      },

      type() {
        if (this.regex === 'scientific') { return 'text'; }
        return 'number';
      },
    },
  };
</script>

<style>
  .settings-label {
    padding-right: 5px;
    width: 50%;
  }

  .settings-input {
    width: 50%;
    height: 2em;
  }

  .settings-row {
    display: flex;
    justify-content: space-between;
    padding-top: 5px;
  }

  .settings-row textarea {
    height: auto;
  }
</style>
