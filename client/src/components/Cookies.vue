<template>
  <div>
    <!-- cookie alert -->
    <b-alert v-model="showCookieBanner" class="position-fixed fixed-bottom m-0 rounded-0"
      style="z-index: 2000;" variant="secondary" dismissible align="center">
      ePotku uses cookies to ensure the best experience during ToF-ERD analysis. <u
        @click="showModal">Click
        here</u> for more information.
    </b-alert>

    <!-- Cookie info modal -->
    <b-modal ref="cookie-info" hide-footer title="ePotku cookie usage">
      <div class="d-block">
        <h3>Some words about cookies</h3>
        <p>
          <ol>
            <li>Users must be informed that cookies are being used on your website, including:</li>
            <ul>
              <li>What cookies are used</li>
              <li>Why they are used</li>
              <li>How they are used</li>
            </ul>
            <li>Prior informed, specific, and voluntary consent must be obtained before cookies are
              placed on a user's computer equipment and before information about a user's computer
              equipment is accessed.When obtaining consent, there are two methods that are allowable
              here:
              <ul>
                <li><strong> Affirmative action/explicit consent</strong></li>

                Clear and explicit affirmative consent can be obtained by placing a check box or a
                clickable button in the notice and requiring a user to click in order to consent.

                The example below shows a button that is labeled with "I Agree" and will work to
                obtain affirmative consent.

                <li><strong>Further browsing/implied consent</strong></li>

                Implied consent will qualify as enough consent to make cookies placement valid so
                long as the following conditions are met:

                There must be a notice that cookies are being used, and it must be displayed in a
                clearly visible and unmissable way on the homepage so that upon first visiting the
                site or using the mobile app, a user will see it
                This notice must make it very clear to the user that by continuing to browse the
                website, consent to place cookies on their device will be implied
                This notice must remain visible until the user actually continues browsing the
                website
                Below is an example of how implied consent can be obtained by using banner ads that
                make it known that continuing to browse will be taken as consent.
              </ul>
            </li>
          </ol>
        </p>
        <p>Words borrowed from <a href="https://www.termsfeed.com/blog/eu-cookies-directive/">here</a> </p>
      </div>
      <b-button class="mt-3" variant="outline-dark" block @click="hideModal">Got It</b-button>
    </b-modal>
  </div>
</template>

<script>
  export default {
    data() {
      return {
        showCookieBanner: true,
      };
    },
    methods: {
      showModal() {
        this.$refs['cookie-info'].show();
      },
      hideModal() {
        this.$refs['cookie-info'].hide();
      },
    },
  };

</script>
