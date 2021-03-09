<template>
  <!-- About page for ePotku -->
  <b-modal id="about" hide-footer title="About ePotku">
    <div class="d-block text-center">
      <h3>ePotku equals e times Potku</h3>
      <b-img
        src="https://www.jyu.fi/science/en/physics/research/infrastructures/accelerator-laboratory/pelletron/potku/banana_logo_small.png/@@images/5d3dd09d-09e0-68fa-ff26-562cb08dd08f.png"
        fluid alt="Responsive image"></b-img>
      <p>
        ePotku is a software project carried out in the Faculty of Information Technologyat the
        University of Jyväskylä in 2020. The project developed a web version of anexisting desktop
        application <a
          href="https://www.jyu.fi/science/en/physics/research/infrastructures/accelerator-laboratory/pelletron/potku/">Potku</a>,
        which is an analysis application developed forthe Accelerator based
        materials physics research group in the Department of Physicsat the University of Jyväskylä.
      </p>
      <p>ePotku was developed by Minja Hänninen, Ilari Jalli, Ville Kuokkanen, Pasi Niininen and
        Tuomas Pitkänen</p>
    </div>
  </b-modal>
</template>

<script>
  export default {
    name: 'About',
  };

</script>
