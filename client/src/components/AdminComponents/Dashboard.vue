<template>
  <!-- Dashboard for admin operations -->
  <div>
    <b-card no-body>
      <b-tabs card active-nav-item-class="text-dark" nav-item-class="bg-dark">
        <!-- Dashboard-tab -->
        <b-tab active>
          <template v-slot:title>
            <b-icon-pie-chart></b-icon-pie-chart>
            Dashboard
          </template>
          <!-- Content to tab -->
          <b-card-text>
            <div class="container-fluid">
              <!-- Heading -->
              <h1 class="h3 mb-2 text-gray-800">Dashboard</h1>

              <p class="mb-4">This is supposed to be dashboard of the Admin page. Here will be
                displayed
                almost-realtime data on usage, logged in users etc.</p>

              <!-- Here are the charts put to an a single row -->
<!--              <div class="row">-->
<!--                <div class="col-xl-8 col-lg-7">-->
<!--                  <div>-->
<!--                    <Plotly class="epotkuUsage" :data="resourceUse" :layout="layoutUsage"-->
<!--                      :display-mode-bar="false" :responsive="false" />-->
<!--                  </div>-->
<!--                </div>-->

                <!-- Pie chart example for biggest users -->
<!--                <div class="col-xl-4 col-lg-5">-->
<!--                  <div class="float">-->
<!--                    <Plotly class="epotkuUsageByUsers"
                          :data="usagebyUsersPie" :layout="layoutUsers"-->
<!--                      :display-mode-bar="false" :responsive="false"/>-->
<!--                  </div>-->
<!--                </div>-->
<!--              </div>-->
            </div>

          </b-card-text>
        </b-tab>

        <!-- Users-tab -->
        <b-tab lazy>
          <template v-slot:title>
            <b-icon-person-fill></b-icon-person-fill>
            Users
          </template>
          <!-- Content to tab -->
          <b-card-text>
            <div class="container-fluid">
              <!-- Heading -->
              <h1 class="h3 mb-2 text-gray-800">Users</h1>
              <Users />
            </div>

          </b-card-text>
        </b-tab>

        <!-- User groups-tab -->
        <b-tab lazy>
          <template v-slot:title>
            <b-icon-people-fill></b-icon-people-fill>
            User Groups
          </template>
          <!-- Content to tab -->
          <b-card-text>
            <div class="container-fluid">
              <!-- Heading -->
              <h1 class="h3 mb-2 text-gray-800">User groups</h1>
              <UserGroup />

            </div>

          </b-card-text>
        </b-tab>

        <!-- Usage tab -->
        <b-tab lazy>
          <template v-slot:title>
            <b-icon-graph-up></b-icon-graph-up>
            Usage
          </template>
          <!-- Content to tab -->
          <b-card-text>
            <div class="container-fluid">
              <!-- Heading -->
              <h1 class="h3 mb-2 text-gray-800">Usage</h1>

              <p class="mb-4">Here's some usage statistics</p>
            </div>

          </b-card-text>
        </b-tab>

        <!-- Log-tab -->
        <b-tab>
          <template v-slot:title>
            <b-icon-archive></b-icon-archive>
            Log
          </template>
          <!-- Content to tab -->
          <b-card-text>
            <div class="container-fluid">
              <!-- Heading -->
              <h1 class="h3 mb-2 text-gray-800">Log</h1>

              <p class="mb-4">Some log information</p>
            </div>
          </b-card-text>
        </b-tab>
      </b-tabs>
    </b-card>
  </div>
</template>
<script>
  // TODO: use plotly.js-cartesian-dist for a smaller build size
  // import { Plotly } from 'vue-plotly';
  import Users from './Users.vue';
  import UserGroup from './UserGroup.vue';

  export default {
    name: 'Dashboard',
    components: {
      Users,
      UserGroup,
      // Plotly,
    },
    data() {
      return {
        resourceUse: [{
          x: ['30.3.', '31.3.', '1.4.', '2.4.', '3.4.', '4.4.', '5.4.'],
          y: [99, 10, 37, 24, 86, 73, 21],
          type: 'scatter',
        }],
        usagebyUsersPie: [{
          type: 'pie',
          labels: ['Seppo Hovi', 'Barack Obama', 'Matti Nykänen', 'Martti Ahtisaari'],
          values: [17, 42, 21, 15],
        }],
        layoutUsers: {
          title: 'Biggest Users',
          showlegend: false,
        },
        layoutUsage: {
          title: 'ePotku Resource Usage %',
        },
      };
    },
  };

</script>
<style>
  .epotkuUsage {
    max-height: 400px;
  }
  .epotkuUsageByUsers {
    max-height: 400px;
  }
</style>
