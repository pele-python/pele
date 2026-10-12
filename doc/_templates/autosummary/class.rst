{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

{% if all_methods %}
.. rubric:: Methods

.. autosummary::
   :toctree:
   {% for item in all_methods %}
      {%- if not item.startswith('_') or item in ['__call__'] %}
   ~{{ objname }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
{% endif %}

{% if all_attributes %}
.. rubric:: Attributes

.. autosummary::
   :toctree:
   {% for item in all_attributes %}
      {%- if not item.startswith('_') %}
   ~{{ objname }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
{% endif %}
