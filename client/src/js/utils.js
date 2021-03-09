/**
 * Shared utility functions
 * TODO: no leading zeroes allowed
 */
export default {
  /**
   * Checks the validity of user inputs and displays a warning if an invalid input is detected
   * @param input Input field whose content is examined
   * @param regex The type of content an input field has
   * @returns {boolean} Returns TRUE if user input was invalid, FALSE otherwise
   */
  checkValidity(input, regex) {
    switch (regex) {
      case 'int':
        if (input.value.match(/.{12,}/)) {
          input.setCustomValidity('Value cannot exceed 11 digits');
          input.reportValidity();
          return true;
        }
        break;
      case 'double':
        if (input.value.match(/([0-9]{12,})|([.][0-9]{4,})/)) {
          input.setCustomValidity('Value cannot exceed 11 significant '
            + 'digits with 3 decimal accuracy');
          input.reportValidity();
          return true;
        }
        break;
      case 'scientific':
        if (input.value.match(/[0-9]{6,}/)) {
          input.setCustomValidity('Values cannot exceed 5 digits');
          input.reportValidity();
          return true;
        }
        if (!input.value.match(/(-?\d+[.]?\d*e-?\d+)/)) {
          input.setCustomValidity('Value must be in the following format: '
            + '3.1e12 or -3.1e-12');
          input.reportValidity();
          return true;
        }
        break;
      default:
        return false;
    }
    return false;
  },

  checkInput(keypress, regex) {
    switch (regex) {
      case 'int':
        // not a number
        if (keypress.key.match(/[^0-9\b]/)) {
          keypress.preventDefault();
        }
        break;
      case 'double':
        // not a number or period
        if (keypress.key.match(/[^0-9\b.]/)) {
          keypress.preventDefault();
        }
        // only one period allowed
        if (keypress.key === '.' && keypress.target.value.match(/[.]/)) {
          keypress.preventDefault();
        }
        if (keypress.target.value.match(/[.](\d{3})/)) {
          keypress.preventDefault();
        }
        break;
      case 'scientific':
        // not e, -, number or period
        if (keypress.key.match(/[^0-9.\be-]/)) {
          keypress.preventDefault();
        }
        // only one period allowed
        if (keypress.key === '.' && keypress.target.value.match(/[.]/)) {
          keypress.preventDefault();
        }
        // only one e allowed
        if (keypress.key === 'e' && keypress.target.value.match(/e/)) {
          keypress.preventDefault();
        }
        // only 2 dashes allowed
        if (keypress.key === '-' && keypress.target.value.match(/-.*-/)) {
          keypress.preventDefault();
        }
        break;
      default: {
        // do nothing
      }
    }
  },
};
