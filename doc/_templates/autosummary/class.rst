{% extends "!autosummary/class.rst" %}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree:
      {% for item in all_methods %}
         {%- if not item.startswith('_') or item in ['__call__'] %}
      ~{{ objname }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree:
      {% for item in all_attributes %}
         {%- if not item.startswith('_') %}
      ~{{ objname }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
   {% endif %}
   {% endblock %}
