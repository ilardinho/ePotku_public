<template>
  <!-- Form for user information -->
  <b-form>
    <b-form-group id="user-id-group">
      <label for="user-id">User ID</label>
      <!-- ID is never editable -->
      <b-form-input id="user-id" :value="user.id" type="number"
        disabled
      />
    </b-form-group>
    <b-form-group id="user-email-group">
      <label for="user-email">Email</label>
      <b-input id="user-email" v-model="user.email" type="email"
        :disabled="isDisabled"
      />
    </b-form-group>
    <!-- TODO: Remove this. Generate a random password on backend and email it instead. -->
    <b-form-group id="user-password-group" v-if="isPasswordEnabled">
      <label for="user-password">Password</label>
      <b-input id="user-password" v-model="user.password" type="password"
        :disabled="isDisabled"
      />
    </b-form-group>
    <b-form-group id="user-first-name-group">
      <label for="user-first-name">First name(s)</label>
      <b-input id="user-first-name" v-model="user.first_name" type="text"
        :disabled="isDisabled"
      />
    </b-form-group>
    <b-form-group id="user-last-name-group">
      <label for="user-last-name">Last name</label>
      <b-input id="user-last-name" v-model="user.last_name" type="text"
        :disabled="isDisabled"
      />
    </b-form-group>
    <b-form-group id="user-is-admin-group">
      <label for="user-is-admin">Is admin</label>
      <b-form-checkbox id="user-is-admin" v-model="user.is_admin" type="boolean"
        :disabled="isDisabled"
      />
    </b-form-group>
  </b-form>
</template>

<script>
  export default {
    props: [
      'user',
      'isDisabled',
      'isPasswordEnabled',
    ],
  };
</script>
