<template>
  <!-- View for the admin page -->
  <div>
    <HeaderNav />
    <div class="container">
      <AdminDashboard />
    </div>
    <Modals />
  </div>
</template>

<script>
  import AdminDashboard from '../components/AdminComponents/Dashboard.vue';
  import HeaderNav from '../components/HeaderNav.vue';
  import Modals from '../components/Modals.vue';

  export default {
    name: 'AdminConsole',
    components: {
      HeaderNav,
      AdminDashboard,
      Modals,
    },
  };

</script>
