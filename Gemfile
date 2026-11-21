source "https://rubygems.org"

# Matches the Jekyll version and plugins GitHub Pages uses.
gem "github-pages", group: :jekyll_plugins

# Needed for `jekyll serve` on Ruby 3+.
gem "webrick"

# Windows / JRuby helpers
platforms :windows, :jruby do
  gem "tzinfo", ">= 1", "< 3"
  gem "tzinfo-data"
end
gem "wdm", "~> 0.1", platforms: [:windows]
