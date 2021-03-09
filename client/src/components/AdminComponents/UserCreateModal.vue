<template>
  <!-- Modal for creating a user -->
  <b-modal id="userCreateModal" title="Create user"
    ok-title="Create user"
    @ok="onSubmit" @cancel="onCancel" @hide="onHide"
  >
    <UserForm id="userForm" :user="this.user" :isDisabled="false" :isPasswordEnabled="true" />
  </b-modal>
</template>

<script>
  import axios from 'axios';

  import UserForm from './UserForm.vue';

  const BACKEND_URL = 'http://localhost:5000/';
  const USER_URL = `${BACKEND_URL}user`;

  export default {
    name: 'UserCreateModal',
    components: {
      UserForm,
    },
    data() {
      return {
        user: {},
      };
    },
    methods: {
      loadForm() {
        this.clearForm();
      },
      clearForm() {
        this.user = {};
      },
      onSubmit() {
        this.postUserData();
        this.clearForm();
      },
      onCancel() {
        this.clearForm();
      },
      onHide() {
        this.clearForm();
      },

      /* HTTP methods */

      // TODO: Display "success" and "failure" messages for HTTP methods

      postUserData() {
        axios.post(USER_URL, this.user, {
          headers: {
            Authorization: this.$store.getters.getAuthHeader,
          },
        })
          .then((response) => {
            // eslint-disable-next-line no-console
            console.log(response);
            // Notify parent about user creation
            this.$emit('update');
          })
          .catch((error) => {
            // eslint-disable-next-line no-console
            console.error(error);
          });
      },
    },
  };
</script>
