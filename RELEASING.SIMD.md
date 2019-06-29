# Main Release

    PIL_BRANCH=5.4.x

    git fetch
    git checkout $PIL_BRANCH
    git branch simd/$PIL_BRANCH

    git checkout simd/split
    git rebase simd/$PIL_BRANCH
    git push -f

    git checkout simd/rgba-convert
    git rebase simd/$PIL_BRANCH
    git push -f

    git checkout simd/resample
    git rebase simd/$PIL_BRANCH
    git push -f

    git checkout simd/filters
    git rebase simd/$PIL_BRANCH
    git push -f

    git checkout simd/box-blur
    git rebase simd/$PIL_BRANCH
    git push -f

    git checkout simd/alpha-composite
    git rebase simd/$PIL_BRANCH
    git push -f

    git checkout simd/color-LUT
    git rebase simd/$PIL_BRANCH
    git push -f

    git checkout simd/info
    git rebase simd/$PIL_BRANCH
    git rm README.rst
    git rebase --continue
    nano ./src/PIL/_version.py  # '5.4.1.post0'
    git rebase --continue
    git push -f

    git checkout simd/$PIL_BRANCH
    git merge simd/split --no-ff
    git merge simd/rgba-convert --no-ff
    git merge simd/resample --no-ff
    git merge simd/filters --no-ff
    git merge simd/box-blur --no-ff
    git merge simd/alpha-composite --no-ff
    git merge simd/info --no-ff

    # Check here

    git tag v$PIL_BRANCH.post0
    ./setup.py sdist --format=gztar
    git push -u uploadcare simd/$PIL_BRANCH
    # set the branch as default on github

# Check Release

  pip uninstall -y pillow; pip uninstall -y pillow-simd; touch ./src/_imaging.c
  CC="ccache cc -msse4" python ./setup.py develop > /dev/null
  nosetests
  pip uninstall -y pillow; pip uninstall -y pillow-simd; touch ./src/_imaging.c
  CC="ccache cc -mavx2" python ./setup.py develop > /dev/null
  nosetests
