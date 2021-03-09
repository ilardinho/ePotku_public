<template>
  <div class="maximize">

    <!-- Overlay to prevent user input while loading/saving is in progress -->
    <b-overlay :show="show" class="maximize">
      <div :id="id" v-on:click.right="undo"/>
    </b-overlay>

    <!-- Selection info modal -->
    <b-modal ref="selection-info" id="selection-info" title="Element selection"
             @ok="saveSelection" size="sm" no-close-on-backdrop>

      <!-- custom buttons on the modal footer -->
      <template v-slot:modal-footer>
        <div class="button-group">
          <b-button variant="danger"    @click="deleteSelection"> Delete </b-button>
          <b-button variant="secondary" @click="cancelSelection"> Cancel </b-button>
          <b-button variant="primary"   @click="saveSelection">   Save   </b-button>
        </div>
      </template>

      <div>
        <label>Type:
          <select v-model="selectionInfo.elementType">
            <option>ERD</option>
            <option>RBS</option>
          </select>
        </label><br/>

        <label>Weight factor:
          <input type="number" step="0.001" v-model="selectionInfo.weightFactor">
        </label><br/>

        <label>Element:
          <input type="input" v-model="selectionInfo.element">
        </label><br/>

        <label>Standard atomic mass:
          <span>mass goes here</span>
        </label><br/>

        <label>Isotope:
          <input v-model="selectionInfo.isotope">
  <!--        <select v-model="selectionInfo.isotope">-->
  <!--          <option>disabled</option>-->
  <!--        </select>-->
        </label><br/>

        <!-- TODO: automatic color selection based on user settings -->
        <label>Selection color:
          <input type="color" id="selection-color" :value="selectionInfo.line.color">
        </label>
      </div>

    </b-modal>

    <!-- Confirmation dialog -->
    <b-modal ref="confirmation-dialog" id="confirmation-dialog" title="Are you sure?"
             @ok="removeAll" size="sm">
      <p>This action cannot be undone.</p>
    </b-modal>

  </div>
</template>

<script>
  import Plotly from 'plotly.js-cartesian-dist';
  import axios from 'axios';
  import robustPointInPolygon from 'robust-point-in-polygon';

  export default {
    name: 'TofeHistogram',
    props: [
      'ticks',
      'selections',
      'tofeData',
      'compressionX',
      'compressionY',
      'xMin',
      'yMin',
      'transpose',
    ],
    data() {
      return {
        id: 'heatmap',
        show: false,
        clickedCoords: [],
        activeSelection: -1,
        okToDraw: false,

        selectionInfo: {
          elementType: 'ERD',
          weightFactor: 1.00,
          element: 'H',
          isotope: 1,
          name: 'New',
          line: { color: '#ff0000' },
          marker: { color: '#ff0000' },
        },

        // ToF-E data and settings below
        // tofeHeatmap: this.traces,
        tofeHeatmap: [{
          z: this.tofeData,
          hoverinfo: 'none',
          showscale: false,
          showlegend: false,
          colorscale: this.colorScale(0.80),
          type: 'heatmap',
          zsmooth: 'best',
          transpose: this.transpose,
        }],
        layout: {
          // autosize: true,
          dragmode: 'pan',
          hoverdistance: 15,
          hovermode: 'closest',
          showlegend: true,
          modebar: {
            color: '#000000',
            activecolor: '#ff0000',
          },
          margin: {
            l: 40,
            r: 0,
            t: 0,
            b: 40,
          },
          legend: {
            itemclick: 'toggleothers',
            xanchor: 'right',
            yanchor: 'top',
            x: 1,
            y: 0.95,
            font: { family: 'Lucida Console', size: 15 },
          },
          sliders: [{
            currentvalue: { visible: false },
            xanchor: 'center',
            yanchor: 'bottom',
            y: 0,
            x: 0.5,
            // pad: {
            //   t: 30,
            // },
            lenmode: 'pixels',
            len: 200,
            active: 6,
            tickcolor: 'transparent',
            font: { color: 'transparent' },
            steps: this.sliderOptions(),
          }],
          xaxis: {
            title: 'Time of Flight (ch)',
            showgrid: false,
            // ticktext is visual only, does not affect coordinates
            ticktext: this.ticks.xText,
            tickvals: this.ticks.xVals,
          },
          yaxis: {
            title: 'Energy (ch)',
            showgrid: false,
            scaleanchor: 'x',
            tickangle: -90,
            // ticktext is visual only, does not affect coordinates
            ticktext: this.ticks.yText,
            tickvals: this.ticks.yVals,
          },
          /* left here in case there's need to use buttons besides the ones in the toolbar */
          updatemenus: [{
            buttons: [{
              label: 'Reload',
              name: 'reload',
              method: 'skip',
            }],
            showactive: false,
            type: 'buttons',
            x: 1,
            xanchor: 'right',
            y: 1,
            yanchor: 'top',
            pad: {
              r: 20,
            },
          }],
          //     args: ['zsmooth', 'best'],
          //     args2: ['zsmooth', false],
          //     label: 'Smoothing',
          //     method: 'restyle',
          //   }],
          //   direction: 'down',
          //   showactive: true,
          //   type: 'buttons',
          //   // pad: { t: 30 },
          //   x: 0,
          //   xanchor: 'left',
          //   y: 1,
          //   yanchor: 'top',
          //   font: { color: '#5072a8' },
          // },
          //   {
          //   buttons: [{
          //     label: 'Selection',
          //     name: 'newSelection',
          //     method: 'skip',
          //   }],
          //   direction: 'down',
          //   pad: { t: 100 },
          //   showactive: true,
          //   active: -1,
          //   type: 'buttons',
          //   x: 0.99,
          //   xanchor: 'right',
          //   y: 1,
          //   yanchor: 'top',
          //   font: { color: '#5072a8' },
          // }, {
          //   buttons: [{
          //     label: 'Undo',
          //     name: 'undo',
          //     method: 'skip',
          //
          //   direction: 'down',
          //   pad: { t: 30 },
          //   showactive: false,
          //   type: 'buttons',
          //   x: 0.99,
          //   xanchor: 'right',
          //   y: 1,
          //   yanchor: 'top',
          //   font: { color: '#5072a8' },
          // }
          // ],
        },
        config: {
          doubleClick: false,
          scrollZoom: true,
          responsive: true,

          // toolbar options below
          displayModeBar: true,
          displaylogo: false,
          // list of button names
          // https://github.com/plotly/plotly.js/blob/master/src/components/modebar/buttons.js
          // list of icon names
          // https://github.com/plotly/plotly.js/blob/master/src/fonts/ploticon.js
          modeBarButtons: [[
            // 'toImage',
            'zoom2d',
            'pan2d',
            {
              name: 'Toggle smoothing',
              icon: Plotly.Icons.plotlylogo,
              click: () => {
                if (this.tofeHeatmap[0].zsmooth === 'best') {
                  this.toggleIconColor('Toggle smoothing', this.layout.modebar.color);
                  Plotly.restyle(this.id, { zsmooth: false }, 0);
                } else {
                  this.toggleIconColor('Toggle smoothing', this.layout.modebar.activecolor);
                  Plotly.restyle(this.id, { zsmooth: 'best' }, 0);
                }
              },
            },
            'resetViews',
            {
              name: 'New selection',
              icon: Plotly.Icons.pencil,
              click: this.toggleSelection,
            }, {
              name: 'Remove all selections',
              icon: Plotly.Icons.eraseshape,
              click: () => {
                this.$refs['confirmation-dialog'].show();
                // this.removeAll,
              },
            }, {
              name: 'Update cuts',
              icon: Plotly.Icons.disk,
              // click: () => { this.$emit('saveCuts'); },
              click: this.saveCuts,
            }]],
        },
      };
    },
    methods: {
      // TODO: separate method for deletions? perhaps communicates the action more clearly
      saveCuts() {
        this.show = true;
        const cuts = [];

        if (this.tofeHeatmap.length > 1) {
          for (let i = 1; i < this.tofeHeatmap.length; i++) {
            // coordinates need to be scaled back up so a shallow copy won't do
            const selection = {
              elementType: this.tofeHeatmap[i].elementType,
              element: this.tofeHeatmap[i].element,
              isotope: this.tofeHeatmap[i].isotope,
              weightFactor: this.tofeHeatmap[i].weightFactor,
              color: this.tofeHeatmap[i].line.color,
              x: [],
              y: [],
              // TODO: scattered element (RBS)
            };
            // skip the closing point since it shouldn't be saved in the backend
            for (let j = 0; j < this.tofeHeatmap[i].x.length - 1; j++) {
              if (this.tofeHeatmap[0].transpose) {
                selection.x.push((this.tofeHeatmap[i].x[j] * this.compressionX) + this.xMin);
                selection.y.push((this.tofeHeatmap[i].y[j] * this.compressionY) + this.yMin);
              } else {
                // flip x/y if transpose == false
                selection.x.push((this.tofeHeatmap[i].y[j] * this.compressionX) + this.xMin);
                selection.y.push((this.tofeHeatmap[i].x[j] * this.compressionY) + this.yMin);
              }
            }
            cuts.push(selection);
          }
        }

        axios.post('http://localhost:5000/cuts', { cuts })
          .then(() => {
            // TODO: display 'cuts saved' message for user
            this.$emit('getCuts');
            this.show = false;
          })
          .catch((error) => {
            // eslint-disable-next-line no-console
            console.error(error);
            this.show = false;
          });
      },


      /**
       * Colorscale for ToF-E heatmap
       * @param multiplier Used for computing the value of each step in the scale
       * @returns {string[][]} Finished colorscale
       */
      colorScale(multiplier) {
        const scale = [
          ['0.0', '#FFFFFF'],
          ['0.1', '#000000'],
          ['0.2', '#000159'],
          ['0.3', '#0005ff'],
          ['0.4', '#00f4ff'],
          ['0.5', '#00ffa6'],
          ['0.6', '#59ff00'],
          ['0.7', '#f6ef00'],
          ['0.8', '#f78c00'],
          ['0.9', '#f70003'],
          ['1.0', '#69001e'],
        ];

        // each step is 0.xx smaller than the one above
        // comment this out to use a linear scale instead
        let step = 1.0;
        for (let i = scale.length - 2; i > 0; i--) {
          step *= multiplier;
          scale[i][0] = step.toString();
        }
        return scale;
      },


      /**
       * Colorscale options on the slider.
       */
      sliderOptions() {
        const options = [{ visible: false }];

        for (let i = 90; i >= 70; i -= 2) {
          options.push({
            label: (i / 100).toString(),
            method: 'restyle',
            args: ['colorscale', [this.colorScale(i / 100)]],
          });
        }

        return options;
      },

      /**
       * Toggle the color of a custom button in the graph toolbar. Also lets caller determine
       * whether traces are allowed to be extended on the graph.
       * Note: in Plotly v1.54.1 there is no way to do this through the API
       */
      toggleIconColor(icon, color, okToDraw) {
        if (okToDraw !== undefined) this.okToDraw = okToDraw;
        const btns = document.getElementsByClassName('modebar-btn');

        for (let i = 0; i < btns.length; i++) {
          if (btns[i].getAttribute('data-title') === icon) {
            // modebar button -> svg -> path -> style -> fill
            btns[i].children[0].children[0].style.fill = color;
          }
        }
      },


      /**
       * Adds an empty trace to the graph and sets activeSelection to the index of the new trace
       * TODO: cancel selection if user switches to another tab?
       */
      toggleSelection() {
        // if a selection is already active, delete current trace and set selection to inactive
        if (this.activeSelection !== -1) {
          this.toggleIconColor('New selection', this.layout.modebar.color, false);
          Plotly.deleteTraces(this.id, [this.activeSelection]);
          this.activeSelection = -1;
        } else {
          this.toggleIconColor('New selection', this.layout.modebar.activecolor, true);
          // default selection info, also used in Application.vue
          const selection = {
            x: [],
            y: [],
            type: 'scatter',
            mode: 'markers+lines',
            hoverinfo: 'none',
            line: { color: '#ff0000' },
            marker: { size: 5, color: '#ff0000' },
            element: 'H',
            name: 'New',
            elementType: 'ERD',
            isotope: 1,
            weightFactor: 1.0,
          };

          Plotly.addTraces((this.id), selection);
          this.activeSelection = this.tofeHeatmap.length - 1;
        }
      },


      /**
       * If a selection is active, undo the last line. If only one or zero points have been
       * placed, also delete the active selection.
       *
       * @param e: mouse right click event
       */
      undo(e) {
        e.preventDefault();
        // check if a selection is active
        if (this.activeSelection > 0) {
          const i = this.activeSelection;

          // if current selection has 1 or 0 points, delete trace
          if (this.tofeHeatmap[i].x.length <= 1) {
            Plotly.deleteTraces(this.id, [i]);
            this.activeSelection = -1;
            this.toggleIconColor('New selection', this.layout.modebar.color, false);
            // changes 'New Selection' button to inactive
            // Plotly.relayout(this.id, { 'updatemenus[1].active': -1 });
            return;
          }

          // remove the last point from trace data
          this.tofeHeatmap[i].x.length -= 1;
          this.tofeHeatmap[i].y.length -= 1;

          Plotly.restyle(
            this.id,
            { x: [this.tofeHeatmap[i].x], y: [this.tofeHeatmap[i].y] },
            this.activeSelection,
          );
        }
      },


      /**
       * Deletes a selection. Not permanent until cuts are saved.
       */
      deleteSelection() {
        Plotly.deleteTraces(this.id, [this.activeSelection]);
        this.$refs['selection-info'].hide();
        this.activeSelection = -1;
        this.toggleIconColor('New selection', this.layout.modebar.color);
      },


      /**
       * Hide the selection modal and do nothing else.
       */
      cancelSelection() {
        this.$refs['selection-info'].hide();
        if (!this.okToDraw) {
          this.activeSelection = -1;
        }
      },


      /**
       * Saves user's input for the selection (element, isotope etc.). Not permanent until
       * cuts are saved.
       */
      saveSelection() {
        // v-model doesn't work with a color picker
        const color = document.getElementById('selection-color').value;

        // update attributes that don't work with v-model
        this.selectionInfo.line.color = color;
        this.selectionInfo.marker.color = color;
        this.selectionInfo.name = this.selectionInfo.isotope.toString().sup()
          + this.selectionInfo.element;

        // this also updates the active selection's attributes to match selectionInfo
        Plotly.update(this.id, this.selectionInfo, this.layout, [this.activeSelection]);

        this.toggleIconColor('New selection', this.layout.modebar.color, false);
        this.activeSelection = -1;
        this.$refs['selection-info'].hide();
      },


      /**
       * Removes all selections from the graph. Not permanent until cuts are saved.
       */
      removeAll() {
        this.tofeHeatmap.length = 1;
        Plotly.update(
          this.id,
          this.tofeHeatmap,
          this.layout,
        );
        this.activeSelection = -1;
        this.toggleIconColor('New selection', this.layout.modebar.color, false);
      },


      /**
       * Updates the selection info template. V-model directly alters the attribute
       * it's bound to so we need a temporary object to apply any changes to in case
       * the user decides to discard their changes.
       */
      copyToTemplate(sel) {
        const { color } = sel.line;
        this.selectionInfo = {
          elementType: sel.elementType,
          weightFactor: sel.weightFactor,
          element: sel.element,
          isotope: sel.isotope,
          line: { color },
          marker: { color },
        };
      },
    },

    // This re-renders the plot every time the user switches to ToF-E tab (if lazy loading is
    // enabled), which has a minor performance impact but avoids a lot of headache elsewhere.
    // For instance, rendering the plot when ToF-E tab is not focused causes the plot to be
    // minimized in the top-left corner. Same thing happens if the browser window is resized while
    // ToF-E tab is inactive.
    // Not lazy loading tabs has a performance impact on the others since the plot is simply hidden
    mounted() {
      for (let i = 0; i < this.selections.length; i++) {
        this.tofeHeatmap.push(this.selections[i]);
      }
      this.tofeHeatmap[0].colorscale = this.colorScale(0.80);

      Plotly.react(
        this.id,
        this.tofeHeatmap,
        this.layout,
        this.config,
      );

      this.toggleIconColor('Toggle smoothing', this.layout.modebar.activecolor);
      const plotDiv = document.getElementById(this.id);

      // click listener for getting coordinates and adding lines to an active selection
      plotDiv.on('plotly_click', (e) => {
        // eslint-disable-next-line prefer-destructuring
        const x = e.points[0].x;
        // eslint-disable-next-line prefer-destructuring
        const y = e.points[0].y;
        // hack to get coordinates for double clicking since plotly_doubleclick does not return data
        this.clickedCoords = [x, y];

        if (e.event.button === 0 && this.okToDraw) { // check for left-click (right click is cancel)
          if (this.activeSelection > 0) {
            // active selection's X/Y coordinate lists
            const activeX = this.tofeHeatmap[this.activeSelection].x;
            const activeY = this.tofeHeatmap[this.activeSelection].y;

            // don't add the same point twice in a row
            if (x === activeX[activeX.length - 1] && y === activeY[activeY.length - 1]) {
              return;
            }

            Plotly.extendTraces(this.id, { x: [[x]], y: [[y]] }, [this.activeSelection]);

            // minimum 3 points before closing a selection
            if (activeX.length >= 3 && x === activeX[0] && y === activeY[0]) {
              this.copyToTemplate(this.tofeHeatmap[this.activeSelection]);
              this.$refs['selection-info'].show();
            }
          }
        }
      });

      // click listener for selections, enables editing on legend double click
      plotDiv.on('plotly_legenddoubleclick', (e) => {
        if (this.activeSelection === -1) {
          this.activeSelection = e.expandedIndex;
          this.copyToTemplate(this.tofeHeatmap[this.activeSelection]);
          this.$refs['selection-info'].show();
        }
        return false;
      });

      // click listener for selections, enables editing when double clicking a selection on the plot
      plotDiv.on('plotly_doubleclick', () => {
        if (this.activeSelection === -1) {
          for (let i = 1; i < this.tofeHeatmap.length; i++) {
            const polygon = [];

            // robustPointInPolygon requires the points to be in [[x0, y0], [x1, y1]] format
            for (let j = 0; j < this.tofeHeatmap[i].x.length; j++) {
              polygon.push([this.tofeHeatmap[i].x[j], this.tofeHeatmap[i].y[j]]);
            }

            // -1 == inside, 0 == boundary, 1 == outside
            if (robustPointInPolygon(polygon, this.clickedCoords) < 1) {
              this.activeSelection = i;
              this.copyToTemplate(this.tofeHeatmap[this.activeSelection]);
              this.$refs['selection-info'].show();
              break;
            }
          }
        }
      });

      /* left here in case there's need to use buttons besides the ones in the toolbar */
      // click listener for buttons
      document.getElementById(this.id).on('plotly_buttonclicked', (menu) => {
        if (menu.button.name === 'reload') { this.$emit('reload'); }
        // if (menu.button.name === 'newSelection') { this.newSelection(); }
        // if (menu.button.name === 'undo') { this.undo(); }
        // if (menu.button.name === 'removeAll') { this.removeAll(); }
      });
    },
    watch: {
      /**
       * Redraw ToF-E graph when its data array changes, i.e. the user has created a new histogram
       * with different compression values
       */
      tofeData() {
        this.tofeHeatmap[0].transpose = this.transpose;
        // if array size has changed, selection coordinate scaling has as well
        this.tofeHeatmap.length = 1;
        Object.keys(this.selections).forEach((key) => {
          this.tofeHeatmap.push(this.selections[key]);
        });

        // data binding reactivity only applies to html elements, so attributes are updated here
        this.tofeHeatmap[0].z = this.tofeData;
        this.layout.xaxis.ticktext = this.ticks.xText;
        this.layout.yaxis.ticktext = this.ticks.yText;
        this.layout.xaxis.tickvals = this.ticks.xVals;
        this.layout.yaxis.tickvals = this.ticks.yVals;

        // without this, the plot doesn't reset to a full view when re-rendered
        if (this.transpose) {
          this.layout.xaxis.range = [0, this.tofeHeatmap[0].z.length];
          this.layout.yaxis.range = [0, this.tofeHeatmap[0].z[0].length];
        } else {
          this.layout.xaxis.range = [0, this.tofeHeatmap[0].z[0].length];
          this.layout.yaxis.range = [0, this.tofeHeatmap[0].z.length];
        }

        Plotly.update(this.id, this.tofeHeatmap, this.layout);
      },
    },
  };
</script>

<style>
  #selection-info label {
    display: flex;
    justify-content: space-between;
  }

  #selection-info select, input, span {
    width: 30%;
  }

  /* TODO: calc() to make this more exact? */
  /* center heatmap toolbar */
  #heatmap .modebar {
    left: 50%;
    transform: translateX(-20%);
  }

  #heatmap {
    width: 100%;
    height: 98%;
    padding-top: 5px;
    padding-bottom: 5px;
  }

  .button-group {
    display: flex;
    margin-left: 100px;
    /*flex-direction: row;*/
    justify-content: space-between;
  }
</style>
