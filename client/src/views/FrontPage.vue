<template>
  <!-- View for the front/ landing page -->
  <div class="FrontPage">
    <HeaderNav />
    <b-row>
      <div class="jumbotron text-center w-100 jumpo">
        <h1>ePotku</h1>
        <p>Potku, analysis software for ToF-ERDA with integrated MCERD</p>
        <p>Create depth profiles from time-of-flight elastic recoil detection analysis
          measurement data. </p>
      </div>
      <div class="container">
        <div class="row ">
          <div class="col-md-4">
            <h3>Introduction – what is Potku?</h3>
            <p>Potku (Finnish: kick, recoil) is an open source analysis software used for creating
              depth profiles from time-of-flight elastic recoil detection analysis measurement data.
              With the Potku 2.0 update, it now also includes MCERD for Monte-Carlo simulations of
              scattered and recoiled particles. The softwarehas a graphical user interface,
              implemented with Python 3.6 (2018), PyQt and Matplotlib. For the analysis itself, a
              backend of standalone command line programs written in C are used. Some of the
              analysis components can also be used independently of the GUI, if necessary.</p>
          </div>
          <div class="col-md-4">
            <h3>Development</h3>
            <p>The development versions are hosted at GitHub at <a
                href="https://github.com/JYU-IBA/potku/tree/development">JYU-IBA/potku</a> (as of
              today 6.8.2019
              Latest release version of Potku 2.0 is in the master branch, but development branch
              might be some bit ahead in few notations).</p>

            <p>The site acts as a hub for development.
              The latest stable release can be downloaded from the downloads section or from the
              GitHub branches. Please see the requirements before installing! If you install Potku
              from GitHub you must compile the C programs using the provided scripts. The files
              available for download also contain the executables.</p>
          </div>
          <div class="col-md-4">
            <h3>Updates</h3>
            <p><b>2020 April (10.4.2020) update!</b> From the beginning of 2020 we have had a IT
              student Juhani Sundell hunting for the multiple bugs and making more permanent fixes
              to the code. Now the isotope efficiency related bug(s) is corrected and the MCERD side
              seems to work quite fluently. Also, the simulated energy spectra are now updating
              instantly when one changes the recoil distribution and the foil-related mismatch
              between slab-analysis and MCERD side is one of the corrected issues.</p>
            <p><a href="https://www.jyu.fi/science/en/physics/research/infrastructures/accelerator-laboratory/pelletron/potku/">Read more...</a></p>
          </div>
        </div>
      </div>
    </b-row>
    <Cookies />
    <Modals />
  </div>
</template>

<script>
  import HeaderNav from '../components/HeaderNav.vue';
  import Modals from '../components/Modals.vue';
  import Cookies from '../components/Cookies.vue';

  export default {
    name: 'FrontPage',
    components: {
      HeaderNav,
      Modals,
      Cookies,
    },
  };

</script>
