# This script removes git tags like 1.2.3a12 both locally and remotely
# https://gist.github.com/matthewmccullough/898798
re_alpha="^[0-9]\.[0-9]\.[0-9]a[0-9]+$"
for tag in `git tag`
do
    if [[ $tag =~  $re_alpha ]]
    then
        echo $tag
        git push origin :$tag
        git tag -d $tag
    fi
done
