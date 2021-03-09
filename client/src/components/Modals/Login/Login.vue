<template>
  <!-- Modal for user login -->
  <!-- TODO: Merge this component with LoginForm? -->
  <b-modal id="login" ref="login-modal" title="Login" ok-title="Login"
    @ok="onSubmit" @submit="onSubmit" @cancel="onCancel" @hide="onHide"
  >
    <LoginForm :credentials="this.credentials" @enterPressed="onEnterPressed" />
  </b-modal>
</template>

<script>
  import LoginForm from './LoginForm.vue';

  export default {
    name: 'Login',
    components: {
      LoginForm,
    },
    data() {
      return {
        credentials: {},
      };
    },
    methods: {
      loadForm() {
        this.clearForm();
      },
      clearForm() {
        this.credentials = {};
      },
      onSubmit() {
        this.$store.dispatch('login', {
          credentials: this.credentials,
        });
        this.clearForm();
      },
      onCancel() {
        this.clearForm();
      },
      onHide() {
        this.clearForm();
      },
      /** Handler for LoginForm's enter event */
      onEnterPressed() {
        this.onSubmit();
        this.$refs['login-modal'].hide();
      },
    },
  };
</script>
