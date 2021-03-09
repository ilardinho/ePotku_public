<template>
  <!-- Form for Login modal -->
  <!-- Handle enter event separately because @submit does not seem to work. -->
  <b-form class="w-100" @keydown.enter="$emit('enterPressed')">
    <b-form-group id="login-email-group">
      <label for="login-email">Email</label>
      <b-input id="login-email" type="email" v-model="credentials.email" />
    </b-form-group>
    <b-form-group id="login-password-group">
      <label for="login-password">Password</label>
      <b-input id="login-password" type="password" v-model="credentials.password" />
    </b-form-group>
  </b-form>
</template>

<script>
  export default {
    props: [
      'credentials',
    ],
  };
</script>
